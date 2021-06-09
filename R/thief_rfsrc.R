#'Temporal hierarchy reconciliation of multivariate Breiman random forests
#'
#'This function fits multivariate Breiman random forests on a multivariate xts timeseries object and then
#'uses ensemble univariate forecast models on all higher levels of temporal aggregation to reconcile the forecasts
#'
#'@importFrom parallel detectCores parLapply makePSOCKcluster setDefaultCluster clusterExport clusterEvalQ stopCluster
#'@importFrom stats ts end start frequency
#'
#'@param y \code{xts matrix}. The outcome series to be modelled. \code{NAs} are currently not supported
#'@param k \code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}. Default
#'is \code{1}, meaning that a final forecast of \code{frequency} horizons will be returned
#'@param lambda \code{numeric}. The Box Cox power transformation parameter for aggregate series. Must be
#'between \code{-1} and \code{2} inclusive. If \code{y_series} contains zeros, \code{lambda} will be set to
#'\code{max(c(0.7, lambda))} to ensure stability of forecasts
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param horizon \code{integer}. The horizon to forecast. Defaults to \code{frequency}
#'@param cores \code{integer}. The number of cores to use
#'@param tune_nodesize \code{logical}. If \code{TRUE}, \code{\link[randomForestSRC]{tune.nodesize}} is used to try and find
#'the optimal nodesize tuning parameter for the multivariate random forest. This can be slow, so the default is \code{FALSE}
#'and a nodesize of \code{10} is used
#'@param predict_quantiles \code{logical}. If \code{TRUE}, a \code{\link[randomForestSRC]{quantreg.rfsrc}} is used
#'to train a second multivariate random forest using quantile loss to predict uncertainty intervals. If \code{FALSE},
#'the distribution of estimates from the \code{2000} original random forests is used to estimate uncertainties
#'@return A \code{list} containing the reconciled forecast distributions for each series in \code{y}. Each element in
#'the \code{list} is a \code{horizon x 1000 matrix} of forecast predictions
#'
#'@seealso \code{\link{ensemble_base}}, \code{\link[randomForestSRC]{rfsrc}},
#'\code{\link[thief]{reconcilethief}}
#'
#'@details Series in \code{y} are aggregated at all possible levels up to annual using \code{\link[thief]{tsaggregates}}.
#'\code{\link[randomForestSRC]{rfsrc}} is used on the unaggregated series, with regressors including exponentially
#'weighted moving averages of each series as well as time features that can be extracted from the series. Run conditions are:
#'\code{ntree = 2000, nsplit = NULL, nodesize = 12}.
#'\code{\link{ensemble_base}} is used on all aggretated levels to find a weighted ensemble of eight
#'univariate forecast models that minimises mean absolute scaled error. Forecasts are then reconciled
#'using \code{\link[thief]{reconcilethief}} and are optionally constrained using non-negative optimisation if there
#'are no negative values in \code{y}. Adjustments to the original unaggregated forecast are incorporated and a
#'distribution of \code{1000} sample paths for each series' forecast are returned
#'
#'@references Athanasopoulos, G., Hyndman, R.,  Kourentzes, N.,  and Petropoulos, F. Forecasting with temporal hierarchies.
#'(2017) European Journal of Operational Research 262(1) 60–74
#'
#'@examples
#'\donttest{
#'library(mvforecast)
#'data("ixodes_vets_dat")
#'
#'#Fit a thief_rfsrc model
#'mod1 <- thief_rfsrc(y = ixodes_vets_dat$y_train,
#'frequency = 52, k = 1,
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
thief_rfsrc = function(y,
                       k = 1,
                       lambda = NULL,
                       frequency = 52,
                       horizon = NULL,
                       predict_quantiles = TRUE,
                       tune_nodesize = FALSE,
                       cores = parallel::detectCores() - 1){

  # Check variables
  if(!xts::is.xts(y)){
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
    # Use rfsrc on the base (unaggregated) series

      cat('\nFitting a rfsrc model to series with frequency', frequency,'\n')

      # Get a dataframe of the series
      rf_data <- data.frame(y)

      # Calculate exponentially weighted moving averages of each series to include as features
      ewma_filter <- function (x, ratio) {
        c(stats::filter(x * ratio, 1 - ratio, "recursive", init = x[1]))
      }

      if(frequency >= 12){
        windows <- unique(ceiling(seq(3, frequency / 2, length.out = 3)))
      } else if (frequency >= 6) {
        windows <- unique(ceiling(seq(1, frequency / 2, length.out = 2)))
      } else {
        windows <- unique(seq(1, frequency, length.out = 2))
      }

      ewmas <- do.call(cbind, lapply(seq_len(ncol(y)), function(x){
        ewma <- matrix(NA, nrow = length(y[,1]), ncol = length(windows))
        for(w in seq_along(windows)){
          ewma[,w] <- suppressWarnings(jitter(ewma_filter(as.vector(zoo::rollmean(y[,x],
                                                                 k = ceiling(windows[w] ^ 0.8),
                                                                 na.pad = TRUE,
                                                                 fill = 0)),
                                         ratio = (2 / (windows[w] + 1))), amount = 0.25))
        }
        ewma
      }))
      rf_data <- cbind(rf_data, ewmas)

      # Train the model
      if(tune_nodesize){
        rf <- randomForestSRC::rfsrc(as.formula(paste0('cbind(',
                                                       paste0(colnames(y), collapse = ','),
                                                       ')~.')),
                                     data = rf_data,
                                     ntree = 1000,
                                     nsplit = NULL,
                                     nodesize = randomForestSRC::tune.nodesize(as.formula(paste0('cbind(',
                                                                                                 paste0(colnames(y), collapse = ','),
                                                                                                 ')~.')),
                                                                               data = rf_data,
                                                                               ntree = 1000,
                                                                               nodesizeTry = c(1:9, seq(10, 100, by = 5)),
                                                                               nsplit = NULL)$nsize.opt)
        opt_nodesize <- rf$nodesize

      } else {
        rf <- randomForestSRC::rfsrc(as.formula(paste0('cbind(',
                                                       paste0(colnames(y), collapse = ','),
                                                       ')~.')),
                                     data = rf_data,
                                     ntree = 1000,
                                     nsplit = NULL,
                                     nodesize = 10)
        opt_nodesize <- 10
      }


      # Store residuals
      for(j in 1:NCOL(y)){
        residuals[[1]][[j]] <- rf$regrOutput[[j]]$predicted - as.vector(y[,j])
      }

      # Generate forecasted ewmas for prediction
      cat('\nForecasting from the rfsrc\n')
      ewma_fcs <- matrix(NA, nrow = frequency * k, ncol = ncol(ewmas))
      for(w in 1:ncol(ewmas)){
        # Add Gaussian noise to forecasted moving averages for better generalizability
        ewma_fcs[,w] <- suppressWarnings(jitter(forecast::snaive(ts(ewmas[,w],
                                                     frequency = frequency),
                                                  h = frequency * k)$mean, amount = 0.25))
      }
      newdata <- data.frame(ewma_fcs)
      colnames(newdata) <- colnames(rf_data)[-c(1:NCOL(y))]

      if(predict_quantiles){
        rf_quantiles <- randomForestSRC::quantreg(as.formula(paste0('cbind(',
                                                                    paste0(colnames(y), collapse = ','),
                                                                    ')~.')),
                                                  data = rf_data,
                                                  ntree = 1000,
                                                  nsplit = NULL,
                                                  nodesize = opt_nodesize,
                                                  method = 'local')
        preds_quantiles <- randomForestSRC::quantreg(object = rf_quantiles , newdata = newdata)
        preds_quantiles <- purrr::map(preds_quantiles$quantreg, 'quantiles')
        rm(rf_quantiles)

      } else {
        # Prediction using the mean forest
        preds_tree <- lapply(1:1000, function(p){
          preds <- randomForestSRC::predict.rfsrc(rf, get.tree = p + 1000, newdata = newdata)
          do.call(cbind, purrr::map(preds$regrOutput, 'predicted'))
        })
      }
      rm(newdata, rf_data, ewmas, ewma_fcs, rf)

      # Loop across series and calculate forecast prediction statistics for later reconciliation
      cat('\nCalculating prediction intervals\n')
      series_forecasts <- lapply(seq_len(ncol(y)), function(series){
        if(predict_quantiles){
          series_quantiles <- preds_quantiles[[series]]
          forecast <- do.call(rbind, lapply(seq_len(frequency * k), function(x){
            quantiles <- c(series_quantiles[x,2], series_quantiles[x,20],
                           series_quantiles[x,50], series_quantiles[x,80],
                           series_quantiles[x,98])
          }))

          # Smooth the uncertainy intervals
          forecast <- do.call(cbind, lapply(seq_len(ncol(forecast)), function(x){
            as.vector(forecast::tsclean(zoo::rollmedian(forecast[,x], k = max(3, floor(horizon / 10)), fill = NA)))
          }))

          # Estimate the distribution at each horizon by assuming it is gaussian and erring on the side
          # of caution to reduce over-confidence
          get_sd = function(median, percentile, value){
            abs(value - median) / abs(qnorm(percentile))
          }

          prediction <- do.call(rbind, lapply(seq_len(nrow(forecast)), function(x){
            dist_sd <- max(unlist(lapply(seq(1:99), function(y){
              ifelse(y!=50, get_sd(series_quantiles[x, 50], y / 100, series_quantiles[x,y]), NA)
            })), na.rm = T)
            # fit <- density(series_quantiles[x,], weights = weights / (sum(weights)))
            rnorm(1000, series_quantiles[x,50], dist_sd)
          }))

          } else {

          prediction <- do.call(cbind,lapply(preds_tree, function(x){
            x[,series]
          }))
          if(!any(y < 0)){
            prediction[prediction < 0] <- 0
          }

          forecast <- do.call(rbind, lapply(seq_len(frequency * k), function(x){
            pred_vals <- prediction[x, ]
            pred_vals <- pred_vals[!is.na(pred_vals)]
            ninetyfives <- suppressWarnings(hpd(pred_vals, 0.95))
            eighties <- suppressWarnings(hpd(pred_vals, 0.8))
            quantiles <- c(ninetyfives[1], eighties[1], mean(prediction[x, ]), eighties[3], ninetyfives[3])
            quantiles
          }))
        }

        # Convert to a forecast class object
        forecast <- data.frame(forecast)
        colnames(forecast) <- c('lower95', 'lower80',
                                'mean', 'upper80', 'upper95')

        lower <- cbind(subset(ts(c(0, forecast$lower80), start = end(tsagg[[1]][[i]]),
                                 frequency = frequencies[i]), start = 2),
                       subset(ts(c(0, forecast$lower95), start = end(tsagg[[1]][[i]]),
                                 frequency = frequencies[i]), start = 2))
        colnames(lower) <- c('80%', '95%')
        mean <- subset(ts(c(0, forecast$mean), start = end(tsagg[[1]][[i]]),
                          frequency = frequencies[i]), start = 2)
        upper <- cbind(subset(ts(c(0, forecast$upper80), start = end(tsagg[[1]][[i]]),
                                 frequency = frequencies[i]), start = 2),
                       subset(ts(c(0, forecast$upper95), start = end(tsagg[[1]][[i]]),
                                 frequency = frequencies[i]), start = 2))
        colnames(upper) <- c('80%', '95%')
        forecast <- list(upper = upper, mean = mean, lower = lower,
                         x = outcomes[[i]][,series], fitted = outcomes[[i]][,series],
                         level = c(80,95),
                         method = 'RF')
        class(forecast) <- 'forecast'
        list(orig_distribution = prediction,
             base_forecast = forecast)
      })

      # Store the original distributions from the base multivariate model for later adjusting
      orig_distributions <- purrr::map(series_forecasts, 'orig_distribution')

      # Extract forecasts into the base list
      for(j in seq_len(ncol(y))){
        base[[1]][[j]] <- .subset2(series_forecasts, j)$base_forecast
      }

    # Use automatic forecasting on aggregated series, in parallel if possible
      if(cores > 1){
        cl <- makePSOCKcluster(cores)
        setDefaultCluster(cl)
        clusterExport(cl, c('frequencies',
                              'outcomes',
                              'lambda',
                              'k',
                              'y'),
                      envir = environment())
        clusterEvalQ(cl, library(forecast))
        clusterEvalQ(cl, library(zoo))
        clusterEvalQ(cl, library(xts))

        cat('\nFitting ensemble forecasts to all remaining series using', cores, 'cores\n')
        ensemble_list <- parLapply(cl, seq(2, length(outcomes)), function(i){
          outcome_base <- vector("list", NCOL(y))
          outcome_residuals <- vector("list", NCOL(y))
          for(j in seq_len(NCOL(y))){

            ensemble <- try(
              suppressWarnings(ensemble_base(y = .subset2(outcomes, i)[,j],
                                             lambda = lambda,
                                             frequency = frequencies[i],
                                             k = k,
                                             bottom_series = FALSE)), silent = TRUE)

            if(inherits(ensemble, 'try-error')){
              rm(ensemble)
              outcome_base[[j]] <- forecast::forecast(.subset2(outcomes, i)[,j],
                                                      h = k * frequencies[i])
              outcome_residuals[[j]] <- residuals(forecast::forecast(.subset2(outcomes, i)[,j],
                                                                     h = k * frequencies[i]))

            } else {
              outcome_base[[j]] <- .subset2(ensemble, 1)
              outcome_residuals[[j]] <- .subset2(ensemble, 2)
              rm(ensemble)
            }

          }
          list(outcome_base = outcome_base, outcome_residuals = outcome_residuals)
        })
        stopCluster(cl)

        for(i in seq(2, length(outcomes))){
          for(j in seq_len(NCOL(y))){
            base[[i]][[j]] <- .subset2(ensemble_list, i-1)$outcome_base[[j]]
          }
        }

        for(i in seq(2, length(outcomes))){
          for(j in seq_len(NCOL(y))){
            residuals[[i]][[j]] <- .subset2(ensemble_list, i-1)$outcome_residuals[[j]]
          }
        }
        rm(ensemble_list)

      } else {

        for(i in seq(2, length(outcomes))){
          cat('\nFitting ensemble forecasts to series at frequency', frequencies[i], '\n')
         for(j in seq_len(NCOL(y))){

           ensemble <- tryCatch({
             suppressWarnings(ensemble_base(y = .subset2(outcomes, i)[,j],
                                            lambda = lambda,
                                            frequency = frequencies[i],
                                            k = k,
                                            bottom_series = ifelse(i == 1, TRUE, FALSE)))
           }, error = function(e) {
             'error'
           })

        if(ensemble == 'error'){
          base[[i]][[j]] <- forecast::forecast(.subset2(outcomes, i)[,j],
                                             h = k * frequencies[i])
          residuals[[i]][[j]] <- residuals(forecast::forecast(.subset2(outcomes, i)[,j],
                                                            h = k * frequencies[i]))

        } else {
          base[[i]][[j]] <- .subset2(ensemble, 1)
          residuals[[i]][[j]] <- .subset2(ensemble, 2)
        }

      }
        }
      }

  # Reconcile the forecasts, use non-negative optimisation constraints if there are no negatives present in y
  cat('\nReconciling original forecasts')
  reconciled <- lapply(seq_len(ncol(y)), function(series){
    series_base <- lapply(seq_along(outcomes), function(x){
      .subset2(base, x)[[series]]
    })
    series_base <- lapply(seq_along(series_base), function(x){
      # In case any forecasts are constant, need to jitter so that covariances can be estimated
      series_base[[x]]$mean <- jitter(series_base[[x]]$mean, amount = 0.001)
      ts(series_base[[x]]$mean, frequency = frequencies[x])
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
    .subset2(series_reconciled, 1)
  })

  # Adjust original distributions using the reconciliation adjustment factors
  adjusted_distributions <- lapply(seq_len(ncol(y)), function(series){
    adjustment <- as.numeric(reconciled[[series]] - ts(base[[1]][[series]]$mean, frequency = frequency))

    new_distribution <- sweep(orig_distributions[[series]], 1, adjustment, "+")
    #new_distribution <- orig_distributions[[series]] + adjustment
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
