#'Temporal hierarchy reconciliation of vector exponential smoothing models
#'
#'This function fits vector exponential smoothing on the base level and optionally on
#'higher levels of temporal aggregation for a multivariate xts timeseries object
#'
#'@importFrom parallel detectCores
#'@importFrom stats ts end start frequency
#'
#'@param y \code{xts matrix}. The outcome series to be modelled. \code{NAs} are currently not supported
#'@param multi_freq \code{integer}. Minimum frequency in the temporal hierarchy that will be modelled using
#'a vector exponential smoothing model. Aggregates with frequencies below this threshold will be modelled with
#'univariate models using the \code{\link[forecast]{forecast}} function
#'@param k \code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}. Default
#'is \code{1}, meaning that a final forecast of \code{frequency} horizons will be returned
#'@param level \code{character}. The dynamics for the level component in the \code{tsvets} model. Options are:
#'c("constant", "diagonal", "common", "full", "grouped")
#'@param slope \code{character}. The dynamics for the slope component in the \code{tsvets} model. Options are:
#'c("none", "constant", "common", "diagonal", "full", "grouped")
#'@param damped \code{character}. The dynamics for the dampening component in the \code{tsvets} model. Options are:
#'c("none", "common", "diagonal", "full", "grouped")
#'@param seasonal \code{character}. The dynamics for the seasonal component in the \code{tsvets} model. Options are:
#' c("none", "common", "diagonal", "full", "grouped")
#'@param group Optional \code{vector} of indices denoting which group the series belongs to (when using the grouped dynamics).
#'Defaults to \code{NULL}
#'@param xreg Optional \code{xts} matrix of external regressors. Defaults to \code{NULL}
#'@param newxreg Optional \code{xts} matrix of future values for external regressors to be used in forecasting.
#'Defaults to \code{NULL}
#'@param xreg_include optional \code{matrix} of dimension \code{ncol(y)} by \code{ncol(xreg)} populated
#'with either 0, 1 or 2+ (0 = no beta, 1 = individual beta and 2 = grouped beta). It is also
#'possible to have group wise pooling. For instance 2 variables sharing one pooled estimates,
#'and 3 other variables sharing another grouped estimate would have values of (2,2,3,3,3).
#'The index for group wise pooling starts at 2 and should be incremented for each new group added. Defaults to \code{NULL}
#'@param lambda \code{numeric proportional}. The multivariate Box Cox power transformation parameter for all series.
#'Must be between \code{0} and \code{1.5} inclusive
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param horizon \code{integer}. The horizon to forecast. Defaults to \code{frequency}
#'@param dependence \code{character}. The multivariate error dependence structure to impose. Options are:
#'c("diagonal", "full", "equicorrelation", "shrinkage")
#'@param cores \code{integer}. The number of cores to use. This is used to initialize the states of each series
#'using \code{\link[tsets]{ets_modelspec}}
#'@param save_plots Logical. Plots of fitted and residual values for the unaggregated series will be saved to
#'\code{fig_path} if \code{TRUE}
#'@param fig_path \code{character}. Optional filepath where fitted and residual plots will be saved. Defaults to the
#'current working directory
#'@return A \code{list} containing the reconciled forecast distributions for each series in \code{y}. Each element in
#'the \code{list} is a \code{horizon x 1000 matrix} of forecast predictions
#'
#'@seealso \code{\link[tsvets]{vets_modelspec}}, \code{\link[tsets]{ets_modelspec}}, \code{\link{ensemble_base}},
#'\code{\link[thief]{reconcilethief}}
#'
#'@details Series in \code{y} are aggregated at all possible levels up to annual using \code{\link[thief]{tsaggregates}}.
#'\code{\link[tsvets]{estimate.tsvets.spec}} is used on the unaggregated (original series), and optionally on higher levels
#'of aggregation down to a frequency of \code{multi_freq}. At frequencies below \code{multi_freq}, the best-fitting
#'univariate model is chosen using automatic ensembles in \code{\link{ensemble_base}}. Forecasts are reconciled
#'using \code{\link[thief]{reconcilethief}} and are optionally constrained using non-negative optimisation if \code{lambda}
#'is provided. Adjustments to the original unaggregated forecast are incorporated and a distribution of \code{1000} sample
#'paths for each series' forecast are returned
#'
#'@references Athanasopoulos, G and de Silva, A. (2012),Multivariate Exponential Smoothing for Forecasting Tourist Arrivals,
#'Journal of Travel Research 51(5) 640–652.\cr\cr
#'de Silva, A., R. Hyndman, and R. D. Snyder. (2010).The Vector Innovations Structural Time Series Framework: A Simple Approach
#'to Multivariate Forecasting, Statistical Modelling (10) 353–74.\cr\cr
#'Athanasopoulos, G., Hyndman, R.,  Kourentzes, N.,  and Petropoulos, F. Forecasting with temporal hierarchies.
#'(2017) European Journal of Operational Research 262(1) 60–74
#'
#'@examples
#'\donttest{
#'library(mvforecast)
#'data("ixodes_vets_dat")
#'
#'# View the returned data
#'head(ixodes_vets_dat$y_train)
#'head(ixodes_vets_dat$xreg_train)
#'ixodes_vets_dat$xreg_include
#'head(ixodes_vets_dat$future_xreg)
#'
#'# Fit a vets model with no regressors and common seasonality with the tsvets package
#'mod1 <- tsvets:::simulate.tsvets.estimate(tsvets:::estimate.tsvets.spec(tsvets:::vets_modelspec(ixodes_vets_dat$y_train,
#'level = "grouped",
#'slope = "none",
#'damped = "none",
#' seasonal = "common",
#' lambda = 1,
#' dependence = "equicorrelation",
#' frequency = 52,
#' cores = parallel::detectCores() - 1,
#' group = ixodes_vets_dat$groups),
#' solver = "solnp",
#' control = list(trace = 0)),
#' nsim = 1000,
#' h = ixodes_vets_dat$h)
#'
#' # Calculate a summary of the summed out-of-sample CRPS
#' calc_crps(simulation = mod1, y_test = ixodes_vets_dat$y_test)
#'
#' # Explore whether reconciliation of temporal hierarchies improves predictions
#' mod2 <- thief_vets(y = ixodes_vets_dat$y_train,
#' multi_freq = 12,
#' level = "grouped",
#' slope = "none",
#' damped = "none",
#' seasonal = "common",
#' lambda = 1,
#' dependence = "equicorrelation",
#' frequency = 52,
#' cores = parallel::detectCores() - 1,
#' group = ixodes_vets_dat$groups,
#' save_plots = FALSE)
#'
#' calc_crps(simulation = mod2, y_test = ixodes_vets_dat$y_test)
#'
#' # Plot one of the forecasts against the true values in the test set
#' plot_vets_preds(simulation = mod2[[4]])
#' points(as.vector(ixodes_vets_dat$y_test[,4]))}
#'
#'@export
#'
thief_vets = function(y,
                      multi_freq = 12,
                      k = 1,
                      level = "grouped",
                      slope = "none",
                      damped = "none",
                      seasonal = "common",
                      lambda = 1,
                      dependence = "equicorrelation",
                      frequency = 52,
                      horizon = NULL,
                      cores = parallel::detectCores() - 1,
                      group = NULL,
                      xreg = NULL,
                      xreg_include = NULL,
                      newxreg = NULL,
                      save_plots = TRUE,
                      fig_path = ''){

  # Check variables
  if (!xts::is.xts(y)) {
    stop("y must be an xts object")
  }

  n <- NCOL(y)
  if (n == 1){
    stop("\ncannot specify a vector ets model with only one series")
  }

  ynames <- colnames(y)
  if(is.null(ynames)) {
    colnames(y) <- paste0("Series", 1:n)
    ynames <- colnames(y)
  }

  level <- match.arg(arg = level[1],
                     choices = c("constant", "diagonal", "common", "full", "grouped"))
  slope <- match.arg(arg = slope[1],
                     choices = c("none", "constant", "common", "diagonal", "full", "grouped"))
  damped <- match.arg(arg = damped[1],
                      choices = c("none", "common", "diagonal", "full", "grouped"))
  seasonal <- match.arg(arg = seasonal[1],
                        choices = c("none", "common", "diagonal", "full", "grouped"))
  dependence <- match.arg(arg = dependence[1],
                          choices = c("diagonal", "full", "equicorrelation", "grouped_equicorrelation", "shrinkage"))

  if (any(c(level, slope, damped, seasonal) %in% "grouped")) {
    if (is.null(group)) stop("\ngroup cannot be NULL for grouped choice")
    if (length(group) != n) stop("\nlength of group vector must be equal to number of cols of y")
    if (max(group) > n) stop("\nmax group > ncol y...check and resubmit")
    if (all(group == max(group))) stop("\ngroup variable is the same for all y. Try common instead")
  } else {
    group <- NULL
  }

  if(!is.null(lambda)){
    if(lambda < 0 || lambda > 1.5) stop('lambda must be between 0 and 1.5 inclusive')
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
  for(i in 1:ncol(y)){
    series <- xts.to.ts(y[, i], freq = frequency)
    series_agg <- thief::tsaggregates(series)
    names <- vector()
    for(j in seq_along(series_agg)){
      names[j] <- paste0('Frequency_', frequency(series_agg[[j]]))
    }
    names(series_agg) <- names
    tsagg[[i]] <- series_agg
  }

  # Put aggregated ys back together in xts format for vets modelling if frequency is
  # greater than or equal to multi_freq. Otherwise, leave as ts format for
  # automatic univariate forecasting
  outcomes <- lapply(seq_along(tsagg[[1]]), function(x){
    freq <- frequency(tsagg[[1]][[x]])

    if(freq >= multi_freq){
      series <- zoo::as.zoo(do.call(cbind, lapply(tsagg, '[[', x)))
      series <- xts::xts(series, lubridate::date_decimal(zoo::index(series)))
      colnames(series) <- colnames(y)
    } else {
      series <- do.call(cbind, lapply(tsagg, '[[', x))
      colnames(series) <- colnames(y)
    }
    series
  })

  # Repeat for xreg if supplied
  if(!is.null(xreg)){
    xregagg <- vector(mode = 'list')
    for(i in 1:ncol(xreg)){
      series <- xts.to.ts(xreg[, i], freq = frequency)
      series_agg <- thief::tsaggregates(series)
      names <- vector()
      for(j in seq_along(series_agg)){
        names[j] <- paste0('Frequency_', frequency(series_agg[[j]]))
      }
      names(series_agg) <- names
      xregagg[[i]] <- series_agg
    }

    xreg_outcomes <- lapply(seq_along(xregagg[[1]]), function(x){
      freq <- frequency(xregagg[[1]][[x]])

      if(freq >= multi_freq){
        series <- zoo::as.zoo(do.call(cbind, lapply(xregagg, '[[', x)))
        series <- xts::xts(series, lubridate::date_decimal(zoo::index(series)))
        colnames(series) <- colnames(xreg)
      } else {
        series <- do.call(cbind, lapply(xregagg, '[[', x))
        colnames(series) <- colnames(xreg)
      }
      series
    })

    # Repeat for future xreg if supplied
    newxregagg <- vector(mode = 'list')
    for(i in 1:ncol(newxreg)){
      series <- xts.to.ts(newxreg[, i], freq = frequency)
      series_agg <- thief::tsaggregates(series)
      names <- vector()
      for(j in seq_along(series_agg)){
        names[j] <- paste0('Frequency_', frequency(series_agg[[j]]))
      }
      names(series_agg) <- names
      newxregagg[[i]] <- series_agg
    }

    newxreg_outcomes <- lapply(seq_along(newxregagg[[1]]), function(x){
      freq <- frequency(newxregagg[[1]][[x]])

      if(freq >= multi_freq){
        series <- zoo::as.zoo(do.call(cbind, lapply(newxregagg, '[[', x)))
        series <- xts::xts(series, lubridate::date_decimal(zoo::index(series)))
        colnames(series) <- colnames(newxreg)
      } else {
        series <- do.call(cbind, lapply(newxregagg, '[[', x))
        colnames(series) <- colnames(newxreg)
      }
      series
    })
  }

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
  vets_frequencies <- which(frequencies >= multi_freq)
  nonvets_frequencies <- which(!frequencies >= multi_freq)
  for(i in vets_frequencies){

      # If xreg supplied and frequency is large enough, fit a VETS model
      if(!is.null(xreg)){
        cat('\nFitting a vets model with regressors to series with frequency', frequencies[i],'\n')
        full_mod <- suppressWarnings(tsvets:::estimate.tsvets.spec(tsvets:::vets_modelspec(outcomes[[i]],
                                                             level = level,
                                                             slope = slope,
                                                             damped = damped,
                                                             seasonal = seasonal,
                                                             lambda = lambda,
                                                             dependence = dependence,
                                                             frequency = frequencies[i],
                                                             cores = cores,
                                                             group = group,
                                                             xreg = xreg_outcomes[[i]],
                                                             xreg_include = xreg_include),
                                              solver = "solnp",
                                              control = list(trace = 0)))

        if(i == 1){
          if(save_plots){

            # Save relevant plots prior to reconciliation
            cat('\nSaving fitted and residual plots in', paste0(fig_path), '\n')
            pdf(paste0(fig_path,'_fitted.pdf'), width = 6, height = 5)
            plot(full_mod)
            dev.off()

            pdf(paste0(fig_path,'_resids.pdf'), width = 6, height = 5)
            plot(full_mod, type = 'residuals')
            dev.off()
          }
        }

        # Generate predictions using supplied future regressors
        p <- tsvets:::predict.tsvets.estimate(full_mod, h = frequencies[i] * k, newxreg = newxreg_outcomes[[i]])

        # If future regressors are not long enough for a forecast of length frequency * k, lengthen the forecast
        p_null <- suppressWarnings(tsvets:::predict.tsvets.estimate(full_mod, h = frequencies[i] * k,
                                           forc_dates = seq.POSIXt(as.POSIXct(end(outcomes[[i]])),
                                                                   length.out = frequencies[i] + 1,
                                                                   by = frequencies[i])[2:(frequencies[i]+1)]))

        # Loop across series and calculate forecast prediction statistics for later reconciliation
        series_forecasts <- lapply(seq_len(ncol(y)), function(series){
          prediction <- t(as.matrix(p$prediction_table$Predicted[[series]]$distribution))
          prediction_null <- t(as.matrix(p_null$prediction_table$Predicted[[series]]$distribution))

          if(nrow(prediction_null) > nrow(prediction)){
            prediction <- rbind(prediction, prediction_null[(nrow(prediction) + 1):nrow(prediction_null),])
          }

          if(!any(y < 0)){
            prediction[prediction < 0] <- 0
          }

          forecast <- do.call(rbind, lapply(seq_len(nrow(prediction)), function(x){
            pred_vals <- prediction[x, ]
            pred_vals <- pred_vals[!is.na(pred_vals)]
            ninetyfives <- suppressWarnings(hpd(pred_vals, 0.95))
            eighties <- suppressWarnings(hpd(pred_vals, 0.8))
            quantiles <- c(ninetyfives[1], eighties[1], mean.default(prediction[x, ]), eighties[3], ninetyfives[3])
            quantiles
          }))

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
                           method = 'VETS')
          class(forecast) <- 'forecast'
          list(orig_distribution = prediction,
               base_forecast = forecast,
               residuals = full_mod$Error[,series])
        })

      } else {
        # If no xreg is supplied, fit a VETS model without regressors
        cat('\nFitting a vets model with no regressors to series with frequency', frequencies[i],'\n')
        full_mod <- suppressWarnings(tsvets:::estimate.tsvets.spec(tsvets:::vets_modelspec(outcomes[[i]],
                                                             level = level,
                                                             slope = slope,
                                                             damped = damped,
                                                             seasonal = seasonal,
                                                             lambda = lambda,
                                                             dependence = dependence,
                                                             frequency = frequencies[i],
                                                             cores = cores,
                                                             group = group),
                                              solver = "solnp",
                                              control = list(trace = 0)))

        if(i == 1){
          if(save_plots){
            # Save relevant plots prior to reconciliation
            cat('\nSaving fitted and residual plots in', paste0(fig_path), '\n')
            pdf(paste0(fig_path,'_fitted.pdf'), width = 6, height = 5)
            plot(full_mod)
            dev.off()

            pdf(paste0(fig_path,'_resids.pdf'), width = 6, height = 5)
            plot(full_mod, type = 'residuals')
            dev.off()
          }
        }

        # Generate forecasts and summarise into forecast objects
        p_null <- suppressWarnings(tsvets:::predict.tsvets.estimate(full_mod, h = frequencies[i] * k,
                                           forc_dates = seq.POSIXt(as.POSIXct(end(outcomes[[i]])),
                                                                   length.out = frequencies[i] + 1,
                                                                   by = frequencies[i])[2:(frequencies[i]+1)]))

        series_forecasts <- lapply(seq_len(ncol(y)), function(series){
          prediction <- t(as.matrix(p_null$prediction_table$Predicted[[series]]$distribution))

          if(!any(y < 0)){
            prediction[prediction < 0] <- 0
          }

          forecast <- do.call(rbind, lapply(seq_len(nrow(prediction)), function(x){
            pred_vals <- prediction[x, ]
            pred_vals <- pred_vals[!is.na(pred_vals)]
            ninetyfives <- suppressWarnings(hpd(pred_vals, 0.95))
            eighties <- suppressWarnings(hpd(pred_vals, 0.8))
            quantiles <- c(ninetyfives[1], eighties[1], ninetyfives[2], eighties[3], ninetyfives[3])
            quantiles
          }))
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
                           method = 'VETS')
          class(forecast) <- 'forecast'
          list(orig_distribution = prediction,
               base_forecast = forecast,
               residuals = full_mod$Error[,series])
        })
      }

      # Store the original distributions from the base multivariate model for later adjusting
      if(i == 1){
        orig_distrubions <- purrr::map(series_forecasts, 'orig_distribution')
      }

      # Extract forecasts and residuals for later reconciliation
      for(j in seq_len(ncol(y))){
        base[[i]][[j]] <- series_forecasts[[j]]$base_forecast
        residuals[[i]][[j]] <- series_forecasts[[j]]$residuals
      }
  }

  if(cores > 1){
    cl <- makePSOCKcluster(cores)
    setDefaultCluster(cl)
    clusterExport(NULL, c('nonvets_frequencies',
                          'outcomes',
                          'lambda',
                          'k',
                          'y'),
                  envir = environment())
    clusterEvalQ(cl, library(forecast))
    clusterEvalQ(cl, library(zoo))
    clusterEvalQ(cl, library(xts))

    cat('\nFitting ensemble forecasts to all remaining series using', cores, 'cores\n')
    ensemble_list <- parLapply(cl, nonvets_frequencies, function(i){
      outcome_base <- list()
      outcome_residuals <- list()
      for(j in seq_len(NCOL(y))){

        ensemble <- try(suppressWarnings(ensemble_base(y_series = outcomes[[i]][,j],
                                                       lambda = lambda,
                                                       y_freq = frequencies[i],
                                                       k = k,
                                                       bottom_series = ifelse(i == 1, TRUE, FALSE))), silent = TRUE)

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
    for(i in nonvets_frequencies){
      for(j in seq_len(NCOL(y))){
        base[[i]][[j]] <- .subset2(ensemble_list, i - max(vets_frequencies))$outcome_base[[j]]
      }
    }

    for(i in nonvets_frequencies){
      for(j in seq_len(NCOL(y))){
        residuals[[i]][[j]] <- .subset2(ensemble_list, i - max(vets_frequencies))$outcome_residuals[[j]]
      }
    }
    rm(ensemble_list)

  } else {

    for(i in nonvets_frequencies){

      # Use automatic forecasting to get best possible result
      cat('\nFitting ensemble forecasts to series at frequency', frequencies[i], '\n')
      for(j in seq_len(NCOL(y))){

        ensemble <- try(suppressWarnings(ensemble_base(y_series = outcomes[[i]][,j],
                                                       lambda = lambda,
                                                       y_freq = frequencies[i],
                                                       k = k,
                                                       bottom_series = ifelse(i == 1, TRUE, FALSE))), silent = TRUE)

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
        series_reconciled <- suppressWarnings(reconcilethief_nonneg(forecasts = series_base,
                                                                        residuals = series_resids,
                                                                        comb = 'struc'))
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
    series_reconciled[[1]]$mean
  })

  # Adjust original distributions using the reconciliation adjustment factors
  adjusted_distributions <- lapply(seq_len(ncol(y)), function(series){
    adjustment <- as.numeric(reconciled[[series]] - base[[1]][[series]]$mean)

    new_distribution <- orig_distrubions[[series]] + adjustment
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
