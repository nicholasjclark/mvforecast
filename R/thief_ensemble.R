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
#'between \code{-1} and \code{2} inclusive
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param horizon \code{integer}. The horizon to forecast. Defaults to \code{frequency}
#'@param cores \code{integer}. The number of cores to use. This is used to initialize the states of each series
#'using \code{\link[tsets]{ets_modelspec}}
#'@param max_agg (optional) \code{integer} specifying the maximum number of temporal aggregation levels
#'to use when reconciling, via the structural scaling method. Useful if higher levels of aggregation
#'are unlikely to have 'seen' recent changes in series dynamics and will likely then result in poor
#'forecasts as a result. Default is \code{NULL}, meaning that all levels of aggregation are used
#'@param discrete \code{logical} Is the series in \code{y} discrete? If \code{TRUE}, use a copula-based method
#'relying on the Probability Integral Transform to map the series to an approximate Gaussian distribution prior to modelling.
#'Forecasts are then back-transformed to the estimated discrete distribution that best fits \code{y}. Default is \code{FALSE}
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
                      cores = parallel::detectCores() - 1,
                      max_agg = NULL,
                      discrete = FALSE){

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

  if(discrete){
    # Store copula details and random draws from each series' estimated discrete distribution
    copula_details <- vector(mode = 'list')
    lambda <- 1
  }

  for(i in 1:NCOL(y)){
    if(NCOL(y) > 1){
      series <- xts.to.ts(y[, i], freq = frequency)
    } else {
      series <- xts.to.ts(y, freq = frequency)
    }

    # Transform to approximate Gaussian if discrete = TRUE
    if(discrete){
      # Convert y to PIT-approximate Gaussian following censoring and NA interpolation
      copula_y <- copula_params(series)

      # The transformed y (approximately Gaussian following PIT transformation)
      series <- copula_y$y_trans

      copula_details[[i]] <- list(copula_y = copula_y,
                                  dist_params = copula_y$params)
    }

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
    clusterEvalQ(cl, library(mvforecast))
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
                                                       discrete = FALSE,
                                                       bottom_series = FALSE)), silent = T)
        if(inherits(ensemble, 'try-error')){
          outcome_base[[j]] <- forecast::forecast(outcomes[[i]][,j],
                                               h = k * frequencies[i])
          outcome_residuals[[j]] <- residuals(forecast::forecast(outcomes[[i]][,j],
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
                                    discrete = FALSE,
                                    bottom_series = FALSE)), silent = F)

          if(inherits(ensemble, 'try-error')){
            base[[i]][[j]] <- forecast::forecast(outcomes[[i]][,j],
                                                 h = k * frequencies[i])
            residuals[[i]][[j]] <- residuals(forecast::forecast(outcomes[[i]][,j],
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
  reconciled <- lapply(seq_len(NCOL(y)), function(series){
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
      orig_resids[is.infinite(orig_resids)] <- NA
      orig_resids <- as.vector(forecast::tsclean(orig_resids))
      orig_resids[is.infinite(orig_resids)] <- NA
      # Resids must be a multiple of frequency for MinT reconciliation
      jitter(tail(orig_resids, floor(length(orig_resids) / frequencies[x]) * frequencies[x]),
             amount = 0.001)
    })

    if(!any(y < 0) & !discrete){
      series_reconciled <- try(suppressWarnings(reconcilethief_restrict(forecasts = series_base,
                                                                  residuals = series_resids,
                                                                  comb = 'sam',
                                                                  max_agg = max_agg,
                                                                  nonnegative = TRUE)),
                               silent = T)
      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- try(suppressWarnings(reconcilethief_restrict(forecasts = series_base,
                                                                        residuals = series_resids,
                                                                        comb = 'struc',
                                                                        max_agg = max_agg,
                                                                        nonnegative = TRUE)),
                                 silent = T)
      }

      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- try(suppressWarnings(reconcilethief_restrict(forecasts = series_base,
                                                                 residuals = series_resids,
                                                                 comb = 'struc')),
                                 silent = T)
      }

    } else {
      series_reconciled <- try(suppressWarnings(reconcilethief_restrict(forecasts = series_base,
                                                               residuals = series_resids,
                                                               max_agg = max_agg,
                                                               comb = 'sam')),
                               silent = T)
      if(inherits(series_reconciled, 'try-error')){
        series_reconciled <- suppressWarnings(reconcilethief_restrict(forecasts = series_base,
                                                             residuals = series_resids,
                                                             max_agg = max_agg,
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
            sd = abs(base[[1]][[series]]$upper[y,2] - base[[1]][[series]]$mean[y]))
    }))
    adjustment <- as.numeric(reconciled[[series]]$mean - base[[1]][[series]]$mean)

    new_distribution <- sweep(orig_distribution, 1, adjustment, "+")

    if(!any(y < 0) & !discrete){
      new_distribution[new_distribution < 0] <- 0
    }

    if(horizon < frequency){
      new_distribution <- new_distribution[1:horizon,]
    }

    if(discrete){
      # Back-transform the predictions to the estimated discrete distribution
      fcast_vec <- as.vector(new_distribution)
      predictions <- back_trans(x = fcast_vec,
                                params = copula_details[[series]]$dist_params)
      out <- matrix(data = predictions, ncol = ncol(new_distribution), nrow = nrow(new_distribution))
    } else {
      out <- new_distribution
    }
    if(any(is.infinite(out))){
      out[is.infinite(out)] <- max(out, na.rm = T)
    }
    out

  })

  # Return the reconciled forecast distributions for each series in y as a list
  names(adjusted_distributions) <- colnames(y)
  return(adjusted_distributions)
}

