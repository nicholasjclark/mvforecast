#'Temporal hierarchy reconciliation of previously fitted univariate forecast
#'
#'This function fits ensemble univariate forecast models on all levels of temporal
#'aggregation for a series and reconciles the supplied univariate forecast object
#'
#'@importFrom stats ts end start frequency
#'
#'@param y Either a \code{xts matrix} or \code{ts} object. The outcome series to be modelled. \code{NAs} are currently not supported
#'@param original_forecast Either a \code{forecast} or \code{list} object, the latter of which must contain slots
#'named \code{forecast} (containing the mean, upper interval and lower interval forecasts) and \code{residuals} (containing
#'residuals for the fitted forecast model)
#'@param forecast_residuals Optional \code{vector} of residuals from the forecast, if the forecast is a matrix representing the
#'forecast distribution
#'@param lambda \code{numeric}. The Box Cox power transformation parameter for all series. Must be
#'between \code{-1} and \code{2} inclusive. If \code{y} contains zeros, \code{lambda} will be set to
#'\code{max(c(0.7, lambda))} to ensure stability of forecasts
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param return_aggregates \code{Logical}. If \code{TRUE}, the forecasts and residuals for the aggregate series will
#'also be returned
#'@param prior_aggregates Optional result of a previous call to \code{thief_reconcile} that contains the
#'forecasts and residuals from the temporal aggregate forecasts for the same \code{y}. If supplied,
#'forecasts for the temporal aggregates will not be calculated again prior to reconciliation, saving on computation.
#'@param max_agg (optional) \code{integer} specifying the maximum number of temporal aggregation levels
#'to use when reconciling, via the structural scaling method. Useful if higher levels of aggregation
#'are unlikely to have 'seen' recent changes in series dynamics and will likely then result in poor
#'forecasts as a result. Default is \code{NULL}, meaning that all levels of aggregation are used
#'
#'@seealso \code{\link[thief]{reconcilethief}}
#'
#'@details Series in \code{y} are aggregated at all possible levels up to annual using \code{\link[thief]{tsaggregates}}.
#'\code{\link{ensemble_base}} is used on all levels of aggregation to find a weighted ensemble of nine
#'univariate forecast models that minimises mean absolute scaled error. Forecasts are then reconciled
#'using \code{\link[thief]{reconcilethief}} and are optionally constrained using non-negative optimisation if there are no
#'negative values in \code{y}. Adjustments to the original supplied forecast are incorporated and
#'this adjusted forecast is returned
#'
#'@references Athanasopoulos, G., Hyndman, R.,  Kourentzes, N.,  and Petropoulos, F. Forecasting with temporal hierarchies.
#'(2017) European Journal of Operational Research 262(1)
#'@examples
#'\donttest{
#'library(mvforecast)
#'data("ixodes_vets_dat")
#'
#'#Fit a univariate model to one of the series in the ixodes data
#'y <- ixodes_vets_dat$y_train[,1]
#'xts.to.ts <- function(x, freq = 52) {
#'start_time <- floor((lubridate::yday(start(x)) / 365) * freq)
#'ts(as.numeric(x), start = c(lubridate::year(start(x)), start_time), freq = freq)}
#'
#'original_forecast <- forecast(auto.arima(xts.to.ts(y, freq = 52)), h = 52)
#'reconciled <- thief_reconcile(y = y, original_forecast = original_forecast, frequency = 52)}
#'
#'#Plot the original and reconciled forecasts
#'reconciled <- thief_reconcile(y = y, original_forecast = original_forecast, frequency = 52)
#'autoplot(original_forecast)
#'autoplot(reconciled)
#'
#'@export
#'
thief_reconcile = function(y,
                           original_forecast,
                           forecast_residuals = NULL,
                           lambda = NULL,
                           frequency,
                           return_aggregates = FALSE,
                           prior_aggregates = NULL,
                           max_agg = NULL){
  if(missing(lambda)){
    lambda <- ifelse((length(which(y == 0)) / length(y)) > 0.1, 0.7, 1)
  }

  lambda <- ifelse(any(as.vector(y) == 0), max(c(0.7, lambda)), lambda)

  if(!is.null(lambda)){
    if(lambda < -1 || lambda > 2) stop('lambda must be between -1 and 2 inclusive')
  }

  xts.to.ts <- function(x, freq = 52) {
    start_time <- floor((lubridate::yday(start(x)) / 365) * freq)
    ts(as.numeric(x),
       start = c(lubridate::year(start(x)),
                 start_time), freq = freq)
  }

  if(class(y) == 'ts'){
    series <- y
  } else {
    series <- xts.to.ts(y, freq = frequency)
  }

  series_agg <- thief::tsaggregates(series)

  names <- vector()
  for(j in seq_along(series_agg)){
    names[j] <- paste0('Frequency_', frequency(series_agg[[j]]))
  }
  names(series_agg) <- names

  # Create objects for storing forecasts and residuals for all aggregates apart from the original
  base <- vector("list", length(series_agg))
  residuals <- vector("list", length(series_agg))

  forecast_class <- class(original_forecast)[1]
  if(forecast_class == 'forecast'){
    base[[1]] <- original_forecast
    if(!(length(original_forecast$mean) / frequency) %in% seq(1,100, by = 1)){
      stop('Forecast horizon must be a multiple of frequency for temporal reconciliation')
    }
    k <- ceiling(length(original_forecast$mean) / frequency)
    residuals[[1]] <- original_forecast$model$residuals
  } else if(forecast_class == 'list'){
    base[[1]] <- original_forecast$forecast
    if(!(length(original_forecast$forecast$mean) / frequency) %in% seq(1,100, by = 1)){
      stop('Forecast horizon must be a multiple of frequency for temporal reconciliation')
    }
    residuals[[1]] <- original_forecast$residuals
    k <- ceiling(length(original_forecast$forecast$mean) / frequency)
  } else {
    quick_fc <- forecast::snaive(y, h = length(apply(original_forecast, 1, mean)))
    quick_fc$mean <- ts(apply(original_forecast, 1, mean), start = start(quick_fc$mean),
                        frequency = frequency)
    base[[1]] <- quick_fc
    if(!(length(apply(original_forecast, 1, mean)) / frequency) %in% seq(1,100, by = 1)){
      stop('Forecast horizon must be a multiple of frequency for temporal reconciliation')
    }
    residuals[[1]] <- forecast_residuals
    quick_fc$residuals <- forecast_residuals
    k <- ceiling(length(apply(original_forecast, 1, mean)) / frequency)
  }

  frequencies <- as.numeric(unlist(lapply(series_agg, frequency), use.names = FALSE))

  if(!is.null(prior_aggregates)){
    for(i in seq(1, length(series_agg) - 1)){
      base[[i+1]] <- prior_aggregates$aggregate_forecasts[[i+1]]
      residuals[[i+1]] <- prior_aggregates$aggregate_residuals[[i+1]]
    }

  } else {

  # Create aggregate forecasts
  for(i in seq(1, length(series_agg) - 1)){
    cat('\nFitting ensemble forecasts to series at frequency', frequencies[i+1], '\n')

    ensemble <- try(suppressWarnings(ensemble_base(y = series_agg[[i+1]],
                                                   lambda = lambda,
                                                   frequency = frequencies[i+1],
                                                   k = k,
                                                   bottom_series = FALSE)), silent = TRUE)

    if(inherits(ensemble, 'try-error')){
      base[[i+1]] <- forecast::forecast(series_agg[[i+1]],
                                        h = k * frequencies[i+1])
      residuals[[i+1]] <- residuals(forecast::forecast(series_agg[[i+1]],
                                                       h = k * frequencies[i+1]))

    } else {
      base[[i+1]] <- ensemble[[1]]
      residuals[[i+1]] <- ensemble[[2]]
    }
  }
  }

  # Reconcile forecasts
  series_base <- lapply(seq_along(base), function(x){
    # In case any forecasts are constant, need to jitter so that covariances can be estimated
    base[[x]]$mean <- jitter(forecast::tsclean(base[[x]]$mean), amount = 0.001)
    base[[x]]
  })
  series_resids <- lapply(seq_along(residuals), function(x){
    orig_resids <- as.vector(forecast::tsclean(residuals[[x]]))
    orig_resids[is.infinite(orig_resids)] <- NA
    # Resids must be a multiple of frequency for MinT reconciliation
    jitter(tail(orig_resids, floor(length(orig_resids) / frequencies[x]) * frequencies[x]),
           amount = 0.001)
  })

  if(!any(y < 0)){
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
  if(forecast_class == 'forecast'){
    reconciled_forecast <- original_forecast
    reconciled_forecast$mean <- series_reconciled[[1]]$mean
    reconciled_forecast$upper <- series_reconciled[[1]]$upper
    reconciled_forecast$lower <- series_reconciled[[1]]$lower
  } else if(forecast_class == 'list'){
    reconciled_forecast <- original_forecast
    reconciled_forecast$forecast$mean <- series_reconciled[[1]]$mean
    reconciled_forecast$forecast$upper <- series_reconciled[[1]]$upper
    reconciled_forecast$forecast$lower <- series_reconciled[[1]]$lower
  } else {
    adjustment <- as.vector(series_reconciled[[1]]$mean - base[[1]]$mean)
    reconciled_forecast <- sweep(original_forecast, 1, adjustment, "+")
    #reconciled_forecast <- original_forecast + adjustment
    if(!any(y < 0)){
      reconciled_forecast[reconciled_forecast < 0] <- 0
    }
  }

  if(return_aggregates){
    output <- list(reconciled_forecast = reconciled_forecast,
                   aggregate_forecasts = series_base,
                   aggregate_residuals = series_resids)
  } else {
    output <- reconciled_forecast
  }

  return(output)
}


