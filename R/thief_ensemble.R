#'Temporal hierarchy reconciliation of ensemble univariate models
#'
#'This function fits ensemble univariate forecast models on all levels of temporal
#'aggregation for a multivariate xts timeseries object
#'
#'@importFrom parallel detectCores parLapply makePSOCKcluster setDefaultCluster clusterExport clusterEvalQ stopCluster
#'@importFrom stats ts end start frequency
#'
#'@param y \code{xts matrix}. The outcome series to be modelled. \code{NAs} are currently not supported
#'@param k \code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}. Default
#'is \code{1}, meaning that a final forecast of \code{frequency} horizons will be returned
#'@param lambda \code{numeric}. The Box Cox power transformation parameter for all series. Must be
#'between \code{-1} and \code{2} inclusive. If \code{y_series} contains zeros, \code{lambda} will be set to
#'\code{max(c(0.7, lambda))} to ensure stability of forecasts
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
#'\code{\link{ensemble_base}} is used on all levels of aggregation to find a weighted ensemble of six
#'univariate forecast models that minimises mean absolute scaled error. Forecasts are then reconciled
#'using \code{\link[thief]{reconcilethief}} and are optionally constrained using non-negative optimisation if there are no
#'negative values in \code{y}. Adjustments to the original unaggregated forecast are incorporated and a distribution of \code{1000} sample
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
#'#Fit a thief_ensemble model
#'mod1 <- thief_ensemble(y = ixodes_vets_dat$y_train,
#'frequency = 52, lambda = 1, k = 1,
#'cores = parallel::detectCores() - 1)
#'
#'#Calculate the out-of-sample CRPS
#'calc_crps(mod1, y_test = ixodes_vets_dat$y_test)
#'
#'Plot simulation results for one of the plots in the NEON dataset
#'plot_mvforecast(simulation = mod1[[4]])
#'points(as.vector(ixodes_vets_dat$y_test[,4]))}
#'
#'@export
#'
thief_ensemble = function(y,
                      k = 1,
                      lambda = NULL,
                      frequency = 52,
                      horizon = NULL,
                      cores = parallel::detectCores() - 1){

  # Check variables
  if (!xts::is.xts(y)) {
    stop("y must be an xts object")
  }

  n <- NCOL(y)
  if(n > 1){
  ynames <- colnames(y)
  if(is.null(ynames)) {
    colnames(y) <- paste0("Series", 1:n)
    ynames <- colnames(y)
  }
  }

  # Automatically set lambda if missing. If many zeros are present in the bottom series, use 0.7, which gives
  # good stability and balances overconfidence of prediction intervals. Otherwise use 1, which does not
  # transform the series but shifts it by -1 (https://otexts.com/fpp2/transformations.html)
  if(missing(lambda)){
    lambda <- ifelse((length(which(y[,1] == 0)) / length(y[,1])) > 0.1, 0.7, 1)
  }

  lambda <- ifelse(any(as.vector(y[,1]) == 0), max(c(0.7, lambda)), lambda)

  if(!is.null(lambda)){
    if(lambda < -1 || lambda > 2) stop('lambda must be between -1 and 2 inclusive')
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
      if(n > 1){
        colnames(series) <- colnames(y)
      } else {
        series <- cbind(series, rep(NA, length(series)))
      }
    series
  })

  # Create objects for storing forecasts and residuals
  base <- vector("list", length(outcomes))
  residuals <- vector("list", length(outcomes))
  for(i in seq_along(outcomes)){
    base[[i]] <- vector("list", NCOL(y))
    residuals[[i]] <- vector("list", NCOL(y))
  }

  # Compute base forecasts using an VETS model if frequency is >= multi_freq, otherwise
  # use automatic forecasting from the forecast package to choose the most appropriate univariate model
  frequencies <- as.numeric(unlist(lapply(tsagg[[1]], frequency), use.names = FALSE))

  if(cores > 1){
    cl <- makePSOCKcluster(cores)
    setDefaultCluster(cl)
    clusterExport(NULL, c('frequencies',
                          'outcomes',
                          'lambda',
                          'k',
                          'y'),
                  envir = environment())
    clusterEvalQ(cl, library(forecast))
    clusterEvalQ(cl, library(zoo))
    clusterEvalQ(cl, library(xts))

    cat('\nFitting ensemble forecasts to all series using', cores, 'cores\n')
    ensemble_list <- parLapply(cl, seq_along(outcomes), function(i){
      outcome_base <- list()
      outcome_residuals <- list()
      for(j in seq_len(NCOL(y))){

        ensemble <- try(suppressWarnings(ensemble_base(y = outcomes[[i]][,j],
                                                       lambda = lambda,
                                                       frequency = frequencies[i],
                                                       k = k,
                                                       bottom_series = ifelse(i == 1, TRUE, FALSE))), silent = TRUE)

        if(inherits(ensemble, 'try-error')){
          outcome_base[[j]] <- forecast::snaive(outcomes[[i]][,j],
                                               h = k * frequencies[i])
          outcome_residuals[[j]] <- residuals(forecast::snaive(outcomes[[i]][,j],
                                                              h = k * frequencies[i]))

        } else {
          outcome_base[[j]] <- ensemble[[1]]
          outcome_residuals[[j]] <- ensemble[[2]]
        }

      }
      list(outcome_base = outcome_base, outcome_residuals = outcome_residuals)
    })
    stopCluster(cl)
    base <- purrr::map(ensemble_list, 'outcome_base')
    residuals <- purrr::map(ensemble_list, 'outcome_residuals')
    rm(ensemble_list)

  } else {

  for(i in seq_along(outcomes)){

      # Use automatic forecasting to get best possible result
      cat('\nFitting ensemble forecasts to series at frequency', frequencies[i], '\n')
      for(j in seq_len(NCOL(y))){

          ensemble <- try(suppressWarnings(ensemble_base(y = outcomes[[i]][,j],
                                    lambda = lambda,
                                    frequency = frequencies[i],
                                    k = k,
                                    bottom_series = ifelse(i == 1, TRUE, FALSE))), silent = TRUE)

          if(inherits(ensemble, 'try-error')){
            base[[i]][[j]] <- forecast::snaive(outcomes[[i]][,j],
                                                 h = k * frequencies[i])
            residuals[[i]][[j]] <- residuals(forecast::snaive(outcomes[[i]][,j],
                                                                h = k * frequencies[i]))

          } else {
            base[[i]][[j]] <- ensemble[[1]]
            residuals[[i]][[j]] <- ensemble[[2]]
          }

      }
  }
  }

  # Reconcile the forecasts, use non-negative optimisation constraints if there are no negatives present in y
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

    if(!any(y < 0)){
      series_reconciled <- try(suppressWarnings(reconcilethief_nonneg(forecasts = series_base,
                                                                  residuals = series_resids,
                                                                  comb = 'sam')),
                               silent = T)
      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- try(suppressWarnings(reconcilethief_nonneg(forecasts = series_base,
                                                                        residuals = series_resids,
                                                                        comb = 'struc')),
                                 silent = T)
      }

      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- try(suppressWarnings(thief::reconcilethief(forecasts = series_base,
                                                                 residuals = series_resids,
                                                                 comb = 'struc')),
                                 silent = T)
      }

    } else {
      series_reconciled <- try(suppressWarnings(thief::reconcilethief(forecasts = series_base,
                                                               residuals = series_resids,
                                                               comb = 'sam')),
                               silent = T)
      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- suppressWarnings(thief::reconcilethief(forecasts = series_base,
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
    if(!any(y < 0)){
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
