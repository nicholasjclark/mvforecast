#'Temporal hierarchy reconciliation of ensemble univariate models
#'
#'This function fits ensemble univariate forecast models on all levels of temporal
#'aggregation for a multivariate xts timeseries object
#'
#'@importFrom parallel detectCores
#'@importFrom stats ts end start frequency
#'
#'@param y \code{xts matrix}. The outcome series to be modelled. \code{NAs} are currently not supported
#'@param k \code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}. Default
#'is \code{1}, meaning that a final forecast of \code{frequency} horizons will be returned
#'@param lambda \code{numeric proportional}. The Box Cox power transformation parameter for all series. Must be
#'between \code{0} and \code{1} inclusive
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param horizon \code{integer}. The horizon to forecast. Defaults to \code{frequency}
#'@param cores \code{integer}. The number of cores to use. This is used to initialize the states of each series
#'using \code{\link[tsets]{ets_modelspec}}
#'@return A \code{list} containing the reconciled forecast distributions for each series in \code{y}. Each element in
#'the \code{list} is a \code{horizon x 1000 matrix} of forecast predictions
#'
#'@seealso \code{\link{ensemble_base}}, \code{\link[forecast]{forecast}},
#'\code{\link[thief]{reconcilethief}}
#'
#'@details Series in \code{y} are aggregated at all possible levels up to annual using \code{\link[thief]{tsaggregates}}.
#'\code{\link{ensemble_base}} is used on all levels of aggregation to find a weighted ensemble of nine
#'univariate forecast models that minimises mean absolute scaled error. Forecasts are then reconciled
#'using \code{\link[thief]{reconcilethief}} and are optionally constrained using non-negative optimisation if \code{lambda}
#'is provided. Adjustments to the original unaggregated forecast are incorporated and a distribution of \code{1000} sample
#'paths for each series' forecast are returned
#'
#'@references Athanasopoulos, G., Hyndman, R.,  Kourentzes, N.,  and Petropoulos, F. Forecasting with temporal hierarchies.
#'(2017) European Journal of Operational Research 262(1) 60–74
#'
#'@examples
#'\donttest{
#'library(mvforecast)
#'data("ixodes_vets_dat")
#'
#'#Fit a a thief_ensemble model
#'mod1 <- thief_ensemble(y = ixodes_vets_dat$y_train,
#'frequency = 52, lambda = 1, k = 1,
#'cores = parallel::detectCores() - 1)
#'
#'#Calculate the out-of-sample CRPS
#'calc_crps(mod1, y_test = ixodes_vets_dat$y_test)
#'
#'Plot simulation results for one of the plots in the NEON dataset
#'plot_vets_preds(simulation = mod1[[4]])
#'points(as.vector(ixodes_vets_dat$y_test[,4]))}
#'
#'@export
#'
thief_ensemble = function(y,
                      k = 1,
                      lambda = 1,
                      frequency = 52,
                      horizon = NULL,
                      cores = parallel::detectCores() - 1){

  # Check variables
  if (!xts::is.xts(y)) {
    stop("y must be an xts object")
  }

  n <- NCOL(y)
  ynames <- colnames(y)
  if(is.null(ynames)) {
    colnames(y) <- paste0("Series", 1:n)
    ynames <- colnames(y)
  }

  if(!is.null(lambda)){
    if(lambda < 0 || lambda > 1) stop('lambda must be between 0 and 1 inclusive')
  }

  # Set forecast horizon if missing
  if(missing(horizon)){
    horizon <- frequency
  }

  # Function to convert xts to ts object
  xts.to.ts <- function(x, freq = 52) {
    start_time <- floor((lubridate::yday(start(x)) / 365) * freq)
    ts(as.numeric(x),
       start = c(lubridate::year(start(x)),
                 start_time), freq = freq)
  }

  # Construct all temporal aggregates for each series in y
  tsagg <- vector(mode = 'list')
  for(i in 1:n){
    series <- xts.to.ts(y[, i], freq = frequency)
    series_agg <- thief::tsaggregates(series)
    names <- vector()
    for(j in seq_along(series_agg)){
      names[j] <- paste0('Frequency_', frequency(series_agg[[j]]))
    }
    names(series_agg) <- names
    tsagg[[i]] <- series_agg
  }

  # Put aggregated ys back together for automatic univariate forecasting
  outcomes <- lapply(seq_along(tsagg[[1]]), function(x){
      series <- do.call(cbind, lapply(tsagg, '[[', x))
      colnames(series) <- colnames(y)
    series
  })

  # Create objects for storing forecasts and residuals
  base <- list()
  residuals <- list()
  for(i in seq_along(outcomes)){
    base[[i]] <- list()
    residuals[[i]] <- list()
  }

  # Compute base forecasts using an VETS model if frequency is >= multi_freq, otherwise
  # use automatic forecasting from the forecast package to choose the most appropriate univariate model
  frequencies <- as.numeric(unlist(lapply(tsagg[[1]], frequency)))
  for(i in seq_along(outcomes)){

      # Use automatic forecasting to get best possible result
      cat('\nFitting ensemble forecasts to series at frequency', frequencies[i], '\n')
      for(j in seq_len(ncol(y))){

          ensemble <- try(suppressWarnings(ensemble_base(y_series = outcomes[[i]][,j],
                                    lambda = lambda,
                                    y_freq = frequencies[i],
                                    k = k,
                                    bottom_series = ifelse(i == 1, TRUE, FALSE))), silent = TRUE)

          if(inherits(ensemble, 'try-error')){
            base[[i]][[j]] <- forecast::forecast(outcomes[[i]][,j],
                                                 lambda = lambda,
                                                 h = k * frequencies[i])
            residuals[[i]][[j]] <- residuals(forecast::forecast(outcomes[[i]][,j],
                                                                lambda = lambda,
                                                                h = k * frequencies[i]))

          } else {
            base[[i]][[j]] <- ensemble[[1]]
            residuals[[i]][[j]] <- ensemble[[2]]
          }

      }
  }

  # Reconcile the forecasts, use non-negative optimisation constraints if lambda is supplied
  cat('\nReconciling original forecasts')
  reconciled <- lapply(seq_len(ncol(y)), function(series){
    series_base <- lapply(seq_along(outcomes), function(x){
      base[[x]][[series]]
    })
    series_base <- lapply(seq_along(series_base), function(x){
      # In case any forecasts are constant, need to jitter so that covariances can be estimated
      series_base[[x]]$mean <- jitter(series_base[[x]]$mean, amount = 0.001)
      series_base[[x]]
    })
    series_resids <- lapply(seq_along(outcomes), function(x){
      orig_resids <- as.vector(residuals[[x]][[series]])
      # Resids must be a multiple of frequency for MinT reconciliation
      jitter(tail(orig_resids, floor(length(orig_resids) / frequencies[x]) * frequencies[x]),
             amount = 0.001)
    })

    if(!is.null(lambda)){
      series_reconciled <- try(suppressWarnings(reconcilethief_nonneg(forecasts = series_base,
                                                                  residuals = series_resids,
                                                                  comb = 'sam')),
                               silent = T)
      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- suppressWarnings(reconcilethief_nonneg(forecasts = series_base,
                                                                        residuals = series_resids,
                                                                        comb = 'struc'))
      }

    } else {
      series_reconciled <- try(suppressWarnings(reconcilethief(forecasts = series_base,
                                                               residuals = series_resids,
                                                               comb = 'sam')),
                               silent = T)
      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- suppressWarnings(reconcilethief(forecasts = series_base,
                                                             residuals = series_resids,
                                                             comb = 'struc'))
      }
    }

    # Return reconciled forecast for the lowest level of aggregation
    list(mean = series_reconciled[[1]]$mean,
         upper = series_reconciled[[1]]$upper,
         lower = series_reconciled[[1]]$lower)
  })

  # Adjust original distributions using the reconciliation adjustment factors
  adjusted_distributions <- lapply(seq_len(ncol(y)), function(series){
    orig_distribution <- do.call(rbind, lapply(seq_len(nrow(base[[1]][[series]]$upper)), function(y){
      rnorm(1000, mean = base[[1]][[series]]$mean[y],
            sd = abs(base[[1]][[series]]$upper[y,2] - base[[1]][[series]]$mean[y]) / 1.46)
    }))
    adjustment <- as.numeric(reconciled[[series]]$mean - base[[1]][[series]]$mean)

    new_distribution <- orig_distribution + adjustment
    if(!is.null(lambda)){
      new_distribution[new_distribution < 0] <- 0
    }

    if(horizon < frequency){
      new_distribution <- new_distribution[1:horizon,]
    }
    new_distribution

  })

  # Return the reconciled forecast distributions for each series in y as a list
  names(adjusted_distributions) <- colnames(y)
  return(adjusted_distributions)
}
