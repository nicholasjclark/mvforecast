#'Fit a series of univariate models and return a weighted ensemble forecast
#'
#'This function fits eight simple univariate forecast models to the series and calculates weights
#'that minimise the mean absolute scaled error of a weighted ensemble forecast
#'
#'@param y_series A \code{ts} object containing the series to forecast
#'@param y_freq The frequency of \code{y_series}
#'@param lambda \code{numeric proportional}. The Box Cox power transformation parameter for all series. Must be
#'between \code{0} and \code{1} inclusive
#'@param k \code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}
#'
#'@details A total of eight simple univariate models are tested on the series. These include
#'\code{\link[forecast]{auto.arima}}, \code{\link[forecast]{tbats}}, \code{\link[forecast]{ets}},
#'\code{\link[forecast]{thetaf}}, \code{\link[forecast]{stlm}} with ARIMA errors, \code{\link[forecast]{rwf}} with drift,
#'\code{\link[forecast]{naive}} and \code{\link[forecast]{snaive}}. The mean absolute scaled error (MASE) of the in-sample
#'series is calculated for each series, and a weighted mean is optimised to minimise in-sample MASE using the
#'"L-BFGS-B" algorithm in \code{\link[stats]{optim}}
#'
#'@return A \code{list} object containing the ensemble forecast and the ensemble residuals
#'
#'@export
ensemble_base = function(y_series, y_freq, lambda, k){

  # If enough observations exist, split into training and testing data
  if(length(y_series) >= (y_freq * 3)){
    cv <- TRUE
    y_train <- subset(y_series, end = floor(length(y_series) * .66))
    y_test <- subset(y_series, start = floor(length(y_series) * .66) + 1)
  } else {
    cv <- FALSE
    y_train <- y_series
    y_test <- rep(NA, (y_freq * k))
  }
  # Try a thetaf model
theta_base <- try(forecast::thetaf(y_train,
                               h = length(y_test)),
                  silent = TRUE)
if(inherits(theta_base, 'try-error')){
  use_theta <- FALSE
  theta_mae <- rep(NA, length(y_train))
} else {
  use_theta <- TRUE
  if(cv){
    theta_mae <- abs(as.vector(theta_base$mean) - as.vector(y_test))
    theta_base <- forecast::thetaf(y_series, h = y_freq * k)
  } else {
    theta_mae <- tail(as.vector(abs(residuals(theta_base))), length(y_test))
  }
}

# Try an auto.arima model
arima_base <- try(forecast::forecast(forecast::auto.arima(y_train,
                                                          lambda = lambda),
                                                          h = length(y_test)),
                                     silent = TRUE)
if(inherits(arima_base, 'try-error')){
  use_arima <- FALSE
  arima_mae <- rep(NA, length(y_train))
} else {
  use_arima <- TRUE
  if(cv){
    arima_mae <- abs(as.vector(arima_base$mean) - as.vector(y_test))
    arima_base <- forecast::forecast(forecast::auto.arima(y_series,
                                                          lambda = lambda),
                                     h = y_freq * k)
  } else {
    arima_mae <- tail(as.vector(abs(residuals(arima_base))), length(y_test))
  }
}

# Try an STLM with arima on the seasonally adjusted series model
stlmar_base <- try(forecast::forecast(forecast::stlm(y_train,
                                                     method = 'arima',
                                                     lambda = lambda),
                                     h = length(y_test)),
                  silent = TRUE)
if(inherits(stlmar_base, 'try-error')){
  use_stlmar <- FALSE
  stlmar_mae <- rep(NA, length(y_train))
} else {
  use_stlmar <- TRUE
  if(cv){
    stlmar_mae <- abs(as.vector(stlmar_base$mean) - as.vector(y_test))
    stlmar_base <- forecast::forecast(forecast::stlm(y_series,
                                                         method = 'arima',
                                                         lambda = lambda),
                                          h = y_freq * k)
  } else {
    stlmar_mae <- tail(as.vector(abs(residuals(stlmar_base))), length(y_test))
  }
}

# Try a TBATS model
tbats_base <- try(forecast::forecast(forecast::tbats(y_train,
                                                 lambda = lambda),
                                   h = length(y_test)),
                silent = TRUE)
if(inherits(tbats_base, 'try-error')){
  use_tbats <- FALSE
  tbats_mae <- rep(NA, length(y_train))
} else {
  use_tbats <- TRUE
  if(cv){
    tbats_mae <- abs(as.vector(tbats_base$mean) - as.vector(y_test))
    tbats_base <- forecast::forecast(forecast::tbats(y_series,
                                                         lambda = lambda),
                                         h = y_freq * k)
  } else {
    tbats_mae <- tail(as.vector(abs(residuals(tbats_base))), length(y_test))
  }
}

# Try an ETS model
if(y_freq <= 24){
  ets_base <- try(forecast::forecast(forecast::ets(y_train,
                                                   lambda = lambda),
                                     h = length(y_test)),
                  silent = TRUE)
} else {
  ets_base <- try(forecast::forecast(forecast::stlf(y_train,
                                                   lambda = lambda),
                                     h = length(y_test)),
                  silent = TRUE)
}

if(inherits(ets_base, 'try-error')){
  use_ets <- FALSE
  ets_mae <- rep(NA, length(y_train))
} else {
  use_ets <- TRUE
  if(cv){
    ets_mae <- abs(as.vector(ets_base$mean) - as.vector(y_test))
    if(y_freq <= 24){
      ets_base <- forecast::forecast(forecast::ets(y_series,
                                                       lambda = lambda),
                                         h = y_freq * k)
    } else {
      ets_base <- forecast::forecast(forecast::stlf(y_series,
                                                   lambda = lambda),
                                     h = y_freq * k)
    }
  } else {
    ets_mae <- tail(as.vector(abs(residuals(ets_base))), length(y_test))
  }
}

# Try a naive model
naive_base <- try(forecast::naive(y_train, lambda = lambda,
                                   h = length(y_test)),
                silent = TRUE)
if(inherits(naive_base, 'try-error')){
  use_naive <- FALSE
  naive_mae <- rep(NA, length(y_train))
} else {
  use_naive <- TRUE
  if(cv){
    naive_mae <- abs(as.vector(naive_base$mean) - as.vector(y_test))
    naive_base <- forecast::naive(y_series, lambda = lambda,
                                      h = y_freq * k)
  } else {
    naive_mae <- tail(as.vector(abs(residuals(naive_base))), length(y_test))
  }
}

# Try an rwf model
rwf_base <- try(forecast::rwf(y_train,
                              lambda = lambda,
                              drift = TRUE,
                              h = length(y_test)),
                  silent = TRUE)
if(inherits(rwf_base, 'try-error')){
  use_rwf <- FALSE
  rwf_mae <- rep(NA, length(y_train))
} else {
  use_rwf <- TRUE
  if(cv){
    rwf_mae <- abs(as.vector(rwf_base$mean) - as.vector(y_test))
    rwf_base <- forecast::rwf(y_series,
                                  lambda = lambda,
                                  drift = TRUE,
                                  h = y_freq * k)
  } else {
    rwf_mae <- tail(as.vector(abs(residuals(rwf_base))), length(y_test))
  }
}

# Try an snaive model
snaive_base <- try(forecast::snaive(y_train,
                                  lambda = lambda,
                                  h = length(y_test)),
                  silent = TRUE)
if(inherits(snaive_base, 'try-error')){
  use_snaive <- FALSE
  snaive_mae <- rep(NA, length(y_train))
} else {
  use_snaive <- TRUE
  if(cv){
    snaive_mae <- abs(as.vector(snaive_base$mean) - as.vector(y_test))
    snaive_base <- forecast::snaive(y_series,
                                        lambda = lambda,
                                        h = y_freq * k)
  } else {
    snaive_mae <- tail(as.vector(abs(residuals(snaive_base))), length(y_test))
  }
}

# Bind MAE estimates together and remove any that are all NAs prior to optimisation
maes <- cbind(theta_mae,
              arima_mae,
              stlmar_mae,
              tbats_mae,
              ets_mae,
              naive_mae,
              rwf_mae,
              snaive_mae)
colnames(maes) <- c('theta', 'arima', 'stlmar', 'tbats', 'ets',
                    'naive', 'rwf', 'snaive')

if(any(which(apply(maes, 2, FUN = function(x) length(which(is.na(x)))) == nrow(maes)))){
  maes <- maes[,-which(apply(maes, 2, FUN = function(x)
    length(which(is.na(x)))) == nrow(maes))]
}

# Scale against the snaive model to convert the MAE to MASE
# The scaling should depend on seasonality and should use naive when there is no seasonality.
# But experimentation suggests the snaive will be
# equivalent to the naive when seasonality is not present
maes <- maes + 1 / (snaive_mae + 1)

# Define the function to minimise
opt_mae <- function(weights, maes) {
  mean(unlist(lapply(seq_len(nrow(maes)), function(h){
    weighted.mean(maes[h,], weights)
  })))
}

# Use optimisation to minimise the function
weights <- rbeta(ncol(maes), 1, 1)
opt <- optim(weights,
             opt_mae,
             maes = maes,
             method = "L-BFGS-B", lower = 0, upper = 1,
             control = list(trace = 0, maxit = 10000))
ens_weights <- opt$par
names(ens_weights) <- colnames(maes)

# Function to return weighted mean forecast object
weight_fcs = function(fcs, ens_weights){
  means <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$mean)}))
  ens_mean <- unlist(lapply(seq_len(nrow(means)), function(x){
    weighted.mean(means[x,], ens_weights)}))

  upper_95s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$upper[,2])}))
  ens_upper95 <- unlist(lapply(seq_len(nrow(upper_95s)), function(x){
    weighted.mean(upper_95s[x,], ens_weights)}))

  upper_80s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$upper[,1])}))
  ens_upper80 <- unlist(lapply(seq_len(nrow(upper_80s)), function(x){
    weighted.mean(upper_80s[x,], ens_weights)}))

  lower_95s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$lower[,2])}))
  ens_lower95 <- unlist(lapply(seq_len(nrow(lower_95s)), function(x){
    weighted.mean(lower_95s[x,], ens_weights)}))

  lower_80s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$lower[,1])}))
  ens_lower80 <- unlist(lapply(seq_len(nrow(lower_80s)), function(x){
    weighted.mean(lower_80s[x,], ens_weights)}))

  ensemble_forecast <- fcs[[1]]
  ensemble_forecast$mean <- ts(ens_mean,
                               start = start(fcs[[1]]$mean), frequency = frequency(fcs[[1]]$mean))
  ensemble_forecast$upper <- cbind(ts(ens_upper80,
                                      start = start(fcs[[1]]$mean), frequency = frequency(fcs[[1]]$mean)),
                                   ts(ens_upper95,
                                      start(fcs[[1]]$mean), frequency = frequency(fcs[[1]]$mean)))
  colnames(ensemble_forecast$upper) <- colnames(fcs[[1]]$upper)

  ensemble_forecast$lower <- cbind(ts(ens_lower80,
                                      start = start(fcs[[1]]$mean), frequency = frequency(fcs[[1]]$mean)),
                                   ts(ens_lower95,
                                      start(fcs[[1]]$mean), frequency = frequency(fcs[[1]]$mean)))
  colnames(ensemble_forecast$lower) <- colnames(fcs[[1]]$lower)
  ensemble_forecast$method <- 'Ensemble'

  # Now get weighted fitted values for calculating residuals
  fitted <- do.call(cbind, lapply(seq_along(fcs), function(x){
      as.vector(fcs[[x]]$fitted)
    }))
  ens_fitted <- unlist(lapply(seq_len(nrow(fitted)), function(x){
    weighted.mean(fitted[x,], ens_weights)}))
  residuals <- ens_fitted - y_series

  return(list(ensemble_forecast, residuals))
}

# Gather the base forecasts and calculate the ensemble
fcs <- list(theta_base,
            arima_base,
            stlmar_base,
            tbats_base,
            ets_base,
            naive_base,
            rwf_base,
            snaive_base)
names(fcs) <- c('theta', 'arima', 'stlmar', 'tbats', 'ets',
                'naive', 'rwf', 'snaive')
fcs <- fcs[names(fcs) %in% names(ens_weights)]
ensemble <- weight_fcs(fcs = fcs, ens_weights = ens_weights)

# Return the ensemble forecast and residuals
return(list(forecast = ensemble[[1]], residuals = ensemble[[2]]))
}
