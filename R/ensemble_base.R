#'Fit a series of univariate models and return a weighted ensemble forecast
#'
#'This function fits nine simple univariate forecast models to the series and calculates weights
#'that minimise the mean absolute scaled error of a weighted ensemble forecast
#'
#'@param y A \code{ts} object containing the series to forecast
#'@param frequency The seasonal frequency of \code{y}
#'@param lambda \code{numeric}. The Box Cox power transformation parameter for all series. Must be
#'between \code{-1} and \code{2} inclusive. If \code{y} contains zeros, \code{lambda} will be set to
#'\code{max(c(0.7, lambda))} to ensure stability of forecasts
#'@param k \code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}
#'@param bottom_series \code{logical}. If \code{TRUE}, the two \code{auto.arima} models will be ignored
#'to ensure bottom level series forecasts are not overconfident
#'@details A total of nine simple univariate models are tested on the series. These include
#'\code{\link[forecast]{auto.arima}} with exponentially weighted moving average regressors,
#'\code{\link[forecast]{auto.arima}} with fourier terms \code{K = 4} regressors,
#'\code{\link[forecast]{tbats}},
#'\code{\link[forecast]{ets}},
#'\code{\link[forecast]{thetaf}},
#'\code{\link[forecast]{stlm}} with ARIMA errors,
#'\code{\link[forecast]{rwf}} with drift,
#'\code{\link[forecast]{naive}} and
#'\code{\link[forecast]{snaive}}.
#'The mean absolute scaled error (MASE) of the in-sample
#'series is calculated for each series, and a weighted mean is optimised to
#'minimise in-sample MASE using the "L-BFGS-B" algorithm in \code{\link[stats]{optim}}.
#'
#'@return A \code{list} object containing the ensemble forecast and the ensemble residuals
#'
#'@export
ensemble_base = function(y, frequency, lambda = NULL, k, bottom_series = FALSE){

  # Automatically set lambda if missing. If many zeros are present, use 0.7, which seems to give
  # good stability and balances overconfidence of prediction intervals. Otherwise use 1, which does not
  # transform the series but shifts it by -1 (https://otexts.com/fpp2/transformations.html)
  if(missing(lambda)){
    lambda <- ifelse((length(which(y == 0)) / length(y)) > 0.1, 0.7, 1)
  }

 lambda <- ifelse(any(as.vector(y) == 0), max(c(0.7, lambda)), lambda)

 if(!is.null(lambda)){
   if(lambda < -1 || lambda > 2) stop('lambda must be between -1 and 2 inclusive')
 }

  cv <- FALSE
  y_train <- y
  y_test <- rep(NA, c(length(y), frequency * k)[which.min(c(0, (frequency * k) - length(y)))])

# Try a thetaf model
theta_base <- try(forecast::thetaf(y_train,
                               h = length(y_test)),
                  silent = TRUE)
if(inherits(theta_base, 'try-error')){
  use_theta <- FALSE
  theta_mae <- rep(NA, length(y_test))
} else {
  use_theta <- TRUE
  if(cv){
    theta_mae <- abs(as.vector(theta_base$mean) - as.vector(y_test))
    theta_base <- forecast::thetaf(y, h = frequency * k)
  } else {
    theta_mae <- tail(as.vector(abs(residuals(theta_base))), length(y_test))
  }
}

if(!bottom_series){

# Try an auto.arima model with exponentially weighted moving average regressors
if(frequency >= 12){
  windows <- unique(ceiling(seq(3, frequency / 2, length.out = 4)))
} else if (frequency >= 6) {
  windows <- unique(ceiling(seq(1, frequency / 2, length.out = 3)))
} else {
  windows <- unique(seq(1, frequency, length.out = 2))
}

ewma_filter <- function (x, ratio) {
  c(stats::filter(x * ratio, 1 - ratio, "recursive", init = x[1]))
}

ewma <- matrix(NA, nrow = length(y_train), ncol = length(windows))
for(i in seq_along(windows)){
  y_unique <- ts(y_train, start = 1, frequency = frequency)

  ewma[,i] <- ewma_filter(as.vector(zoo::rollmean(y_unique,
                                                  k = ceiling(windows[i] ^ 0.8),
                                                  na.pad = TRUE,
                                                  fill = 0)),
                          ratio = (2 / (windows[i] + 1)))
}

ewma_fc <- matrix(NA, nrow = frequency * k, ncol = ncol(ewma))
for(i in 1:ncol(ewma)){
  # Add Gaussian noise to forecasted moving averages for better generalizability
  ewma_fc[,i] <- jitter(forecast::snaive(ts(ewma[,i],
                                              frequency = frequency), h = frequency * k)$mean, amount = 0.25)
}

arima_base <- try(forecast::forecast(forecast::auto.arima(y_train,
                                                          lambda = lambda,
                                                          xreg = ewma),
                                     h = length(y_test),
                                     xreg = ewma_fc),
                                     silent = TRUE)

if(inherits(arima_base, 'try-error')){
  use_arima <- FALSE
  arima_mae <- rep(NA, length(y_test))
} else {
  use_arima <- TRUE
  if(cv){
    arima_mae <- abs(as.vector(arima_base$mean) - as.vector(y_test))

    ewma <- matrix(NA, nrow = length(y_train), ncol = length(windows))
    for(i in seq_along(windows)){
      ewma[,i] <- ewma_filter(as.vector(zoo::rollmean(y_train,
                                                      k = ceiling(windows[i] ^ 0.8),
                                                      na.pad = TRUE,
                                                      fill = 0)),
                              ratio = (2 / (windows[i] + 1)))
    }

    ewma_fc <- matrix(NA, nrow = frequency * k, ncol = ncol(ewma))
    for(i in 1:ncol(ewma)){
      ewma_fc[,i] <- jitter(forecast::snaive(ewma[,i], h = frequency * k)$mean, amount = 0.1)
    }

    arima_base <- forecast::forecast(forecast::auto.arima(y_train,
                                                          lambda = lambda,
                                                          xreg = ewma),
                                     h = frequency * k,
                                     xreg = ewma_fc)
  } else {
    arima_mae <- tail(as.vector(abs(residuals(arima_base))), length(y_test))
  }
}

# Try an auto.arima model with four harmonic Fourier terms as regressors
arimaf_base <- try(forecast::forecast(forecast::auto.arima(y_train,
                                                           lambda = lambda,
                                                           xreg = forecast::fourier(y_train, K = 4)),
                                     h = length(y_test),

                                     # Add Gaussian noise to forecasted terms for better generalizability
                                     xreg = jitter(forecast::fourier(y_train, K = 4, h = length(y_test)),
                                                   amount = 0.1)),
                  silent = TRUE)
if(inherits(arimaf_base, 'try-error')){
  use_arimaf <- FALSE
  arimaf_mae <- rep(NA, length(y_test))
} else {
  use_arimaf <- TRUE
  if(cv){
    arimaf_mae <- abs(as.vector(arimaf_base$mean) - as.vector(y_test))
    arimaf_base <- forecast::forecast(forecast::auto.arima(y,
                                                           lambda = lambda,
                                                           xreg = forecast::fourier(y, K = 4)),
                                      h = frequency * k,
                                      xreg = jitter(forecast::fourier(y, K = 4, h = frequency * k),
                                                    amount = 0.1))
  } else {
    arimaf_mae <- tail(as.vector(abs(residuals(arimaf_base))), length(y_test))
  }
}
}

# Try an STLM with arima on the seasonally adjusted series model
stlmar_base <- try(forecast::forecast(forecast::stlm(y_train,
                                                     method = 'arima',
                                                     lambda = lambda,
                                                     s.window = 'periodic'),
                                     h = length(y_test)),
                  silent = TRUE)
if(inherits(stlmar_base, 'try-error')){
  use_stlmar <- FALSE
  stlmar_mae <- rep(NA, length(y_test))
} else {
  use_stlmar <- TRUE
  if(cv){
    stlmar_mae <- abs(as.vector(stlmar_base$mean) - as.vector(y_test))
    stlmar_base <- forecast::forecast(forecast::stlm(y,
                                                         method = 'arima',
                                                     s.window = 'periodic'),
                                          h = frequency * k)
  } else {
    stlmar_mae <- tail(as.vector(abs(residuals(stlmar_base))), length(y_test))
  }
}

# Try a TBATS model
tbats_base <- try(forecast::forecast(forecast::tbats(y_train,
                                                 lambda = lambda),
                                   h = length(y_test),
                                   lambda = lambda),
                silent = TRUE)
if(inherits(tbats_base, 'try-error')){
  use_tbats <- FALSE
  tbats_mae <- rep(NA, length(y_test))
} else {
  use_tbats <- TRUE
  if(cv){
    tbats_mae <- abs(as.vector(tbats_base$mean) - as.vector(y_test))
    tbats_base <- forecast::forecast(forecast::tbats(y,
                                                         lambda = lambda),
                                         h = frequency * k,
                                     lambda = lambda)
  } else {
    tbats_mae <- tail(as.vector(abs(residuals(tbats_base))), length(y_test))
  }
}

# Try an ETS model
if(frequency <= 24){
  ets_base <- try(forecast::forecast(forecast::ets(y_train, lambda = lambda),
                       h = length(y_test),
                       lambda = lambda), silent = T)

} else {
  ets_base <- try(forecast::forecast(forecast::stlf(y_train,
                                      lambda = lambda),
                       h = length(y_test),
                       lambda = lambda), silent = T)
}

if(inherits(ets_base, 'try-error')){
  use_ets <- FALSE
  ets_mae <- rep(NA, length(y_test))
} else {
  use_ets <- TRUE
  if(cv){
    ets_mae <- abs(as.vector(ets_base$mean) - as.vector(y_test))
    if(frequency <= 24){
      ets_base <- forecast::forecast(ets(y,
                                         h = frequency * k,
                                         lambda = lambda))
    } else {
      ets_base <- forecast::forecast(forecast::stlf(y,
                                                   lambda = lambda),
                                     h = frequency * k,
                                     lambda = lambda)
      ets_base$upper <- cbind(ets_base$upper, ets_base$upper)
      ets_base$lower <- cbind(ets_base$lower, ets_base$lower)
    }
  } else {
    ets_mae <- tail(as.vector(abs(residuals(ets_base))), length(y_test))
  }
}

# Try a naive model
naive_base <- try(forecast::naive(y_train,
                                   h = length(y_test)),
                silent = TRUE)
if(inherits(naive_base, 'try-error')){
  use_naive <- FALSE
  naive_mae <- rep(NA, length(y_test))
} else {
  use_naive <- TRUE
  if(cv){
    naive_mae <- abs(as.vector(naive_base$mean) - as.vector(y_test))
    naive_base <- forecast::naive(y,
                                      h = frequency * k)
  } else {
    naive_mae <- tail(as.vector(abs(residuals(naive_base))), length(y_test))
  }
}

# Try an rwf model
rwf_base <- try(forecast::rwf(y_train,
                              drift = TRUE,
                              h = length(y_test)),
                  silent = TRUE)
if(inherits(rwf_base, 'try-error')){
  use_rwf <- FALSE
  rwf_mae <- rep(NA, length(y_test))
} else {
  use_rwf <- TRUE
  if(cv){
    rwf_mae <- abs(as.vector(rwf_base$mean) - as.vector(y_test))
    rwf_base <- forecast::rwf(y,
                                  drift = TRUE,
                                  h = frequency * k)
  } else {
    rwf_mae <- tail(as.vector(abs(residuals(rwf_base))), length(y_test))
  }
}

# Try an snaive model
snaive_base <- try(forecast::snaive(y_train,
                                  h = length(y_test)),
                  silent = TRUE)
if(inherits(snaive_base, 'try-error')){
  use_snaive <- FALSE
  snaive_mae <- rep(NA, length(y_test))
} else {
  use_snaive <- TRUE
  if(cv){
    snaive_mae <- abs(as.vector(snaive_base$mean) - as.vector(y_test))
    snaive_base <- forecast::snaive(y,
                                        h = frequency * k)
  } else {
    snaive_mae <- tail(as.vector(abs(residuals(snaive_base))), length(y_test))
  }
}

# Bind MAE estimates together and remove any that are all NAs prior to optimisation
if(!bottom_series){
  maes <- cbind(theta_mae,
                arima_mae,
                arimaf_mae,
                stlmar_mae,
                tbats_mae,
                ets_mae,
                naive_mae,
                rwf_mae,
                snaive_mae)
  colnames(maes) <- c('theta', 'arima', 'arimaf', 'stlmar', 'tbats', 'ets',
                      'naive', 'rwf', 'snaive')
} else {
  maes <- cbind(theta_mae,
                stlmar_mae,
                tbats_mae,
                ets_mae,
                naive_mae,
                rwf_mae,
                snaive_mae)
  colnames(maes) <- c('theta', 'stlmar', 'tbats', 'ets',
                      'naive', 'rwf', 'snaive')
}


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
  mean.default(unlist(lapply(seq_len(nrow(maes)), function(h){
    weighted.mean(maes[h,], weights)
  }), use.names = FALSE))
}

# Use optimisation to minimise the function
weights <- rbeta(ncol(maes), 1, 1)
opt <- optim(weights,
             opt_mae,
             maes = maes,
             method = "L-BFGS-B", lower = 0, upper = 1,
             control = list(trace = 0, maxit = 100000))
ens_weights <- opt$par
names(ens_weights) <- colnames(maes)

# Function to return weighted mean forecast object
weight_fcs = function(fcs, ens_weights){
  means <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$mean)}))
  ens_mean <- unlist(lapply(seq_len(nrow(means)), function(x){
    weighted.mean(means[x,], ens_weights)}), use.names = F)
  rm(means)

  upper_95s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$upper[,2])}))
  ens_upper95 <- unlist(lapply(seq_len(nrow(upper_95s)), function(x){
    weighted.mean(upper_95s[x,], ens_weights)}), use.names = F)
  rm(upper_95s)

  upper_80s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$upper[,1])}))
  ens_upper80 <- unlist(lapply(seq_len(nrow(upper_80s)), function(x){
    weighted.mean(upper_80s[x,], ens_weights)}), use.names = F)
  rm(upper_80s)

  lower_95s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$lower[,2])}))
  ens_lower95 <- unlist(lapply(seq_len(nrow(lower_95s)), function(x){
    weighted.mean(lower_95s[x,], ens_weights)}), use.names = F)
  rm(lower_95s)

  lower_80s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$lower[,1])}))
  ens_lower80 <- unlist(lapply(seq_len(nrow(lower_80s)), function(x){
    weighted.mean(lower_80s[x,], ens_weights)}), use.names = F)
  rm(lower_80s)

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
  rm(fcs)
  ens_fitted <- unlist(lapply(seq_len(nrow(fitted)), function(x){
    weighted.mean(fitted[x,], ens_weights)}), use.names = F)
  residuals <- ens_fitted - y
  residuals[is.na(residuals)] <- 0

  return(list(ensemble_forecast, residuals))
}

# Gather the base forecasts and calculate the ensemble
if(!bottom_series){
  fcs <- list(theta_base,
              arima_base,
              arimaf_base,
              stlmar_base,
              tbats_base,
              ets_base,
              naive_base,
              rwf_base,
              snaive_base)
  names(fcs) <- c('theta', 'arima', 'arimaf','stlmar', 'tbats', 'ets',
                  'naive', 'rwf', 'snaive')
} else {
  fcs <- list(theta_base,
              stlmar_base,
              tbats_base,
              ets_base,
              naive_base,
              rwf_base,
              snaive_base)
  names(fcs) <- c('theta', 'stlmar', 'tbats', 'ets',
                  'naive', 'rwf', 'snaive')
}

fcs <- fcs[names(fcs) %in% names(ens_weights)]
ensemble <- weight_fcs(fcs = fcs, ens_weights = ens_weights)
rm(fcs, ens_weights, maes,
   theta_base,
   stlmar_base,
   tbats_base,
   ets_base,
   naive_base,
   rwf_base,
   snaive_base)

# Return the ensemble forecast and residuals
return(list(forecast = .subset2(ensemble, 1), residuals = .subset2(ensemble, 2)))
}
