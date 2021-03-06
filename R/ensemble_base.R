#'Fit a series of univariate models and return a weighted ensemble forecast
#'
#'This function fits six simple univariate forecast models to the series and calculates weights
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
#'@param discrete \code{logical} Is the series in \code{y} discrete? If \code{TRUE}, use a copula-based method
#'relying on the Probability Integral Transform to map the series to an approximate Gaussian distribution prior to modelling.
#'Forecasts are then back-transformed to the estimated discrete distribution that best fits \code{y}. Default is \code{FALSE}
#'@details A total of seven simple univariate models are tested on the series. These include
#'\code{\link[forecast]{stlf}} with an AR model for the seasonally-adjusted series,
#'\code{\link[forecast]{ets}},
#'\code{\link[forecast]{auto.arima}},
#'\code{\link[forecast]{auto.arima}} with fourier terms \code{K = 4} regressors,
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
ensemble_base = function(y, frequency, lambda = NULL, k = 1, bottom_series = FALSE, discrete = FALSE){

  # Check variables
  if (!is.ts(y)) {
    stop("y must be an ts object")
  }

  # Transform to approximate Gaussian if discrete = TRUE
  if(discrete){
    # Convert y to PIT-approximate Gaussian following censoring and NA interpolation
    copula_y <- copula_params(y, non_neg = T, censor = 0.99)

    # Estimate copula parameters from most recent values of y so the
    # returned discrete distribution is more reflective of recent history
    dist_params <- copula_params(tail(y, min(length(y), 100)),
                                 non_neg = T, censor = 0.99)$params

    # The transformed y (approximately Gaussian following PIT transformation)
    y <- copula_y$y_trans

    # Set BoxCox transformation parameter to 1 for no transformation
    lambda <- 1

    # Take random draws from the estimated discrete distribution so predictions can be mapped
    # back to the original data distribution
    if(length(dist_params) == 2){
      dist_mappings <- stats::rnbinom(5000, size = dist_params[1],
                                      mu = dist_params[2])
    } else {
      dist_mappings <- stats::rpois(5000, lambda = dist_params)
    }
  }

  # Automatically set lambda if missing
  if(missing(lambda)){
    if(any(y < 0)){
      lambda <- 1
    } else {
      lambda <- forecast::BoxCox.lambda(y)
    }
    lambda <- ifelse(any(as.vector(y) == 0), max(c(0.7, lambda)), lambda)
  }

 if(!is.null(lambda)){
   if(lambda < -1 || lambda > 2) stop('lambda must be between -1 and 2 inclusive')
 }

  cv <- FALSE
  y_train <- y
  y_test <- rep(NA, c(length(y), frequency * k)[which.min(c(0, (frequency * k) - length(y)))])

  # Try a stlf model
  stlf_base <- try(forecast::forecast(forecast::stlm(y_train, modelfunction = ar),
                                      h = length(y_test)),
                   silent = TRUE)
  if(inherits(stlf_base, 'try-error')){
    use_stlf <- FALSE
    stlf_mae <- rep(NA, length(y_test))
  } else {
    use_stlf <- TRUE
    if(cv){
      stlf_mae <- abs(as.vector(stlf_base$mean) - as.vector(y_test))
      stlf_base <- forecast::forecast(forecast::stlm(y, modelfunction = ar),
                                      h = frequency * k)
    } else {
      stlf_mae <- tail(as.vector(abs(residuals(stlf_base))), length(y_test))
    }
  }

  # Try an ets model
  if(frequency <= 24){

    ets_base <- try(forecast::forecast(forecast::ets(y_train),
                                       h = length(y_test)),
                    silent = TRUE)
    if(inherits(ets_base, 'try-error')){
      use_ets <- FALSE
      ets_mae <- rep(NA, length(y_test))
    } else {
      use_ets <- TRUE
      if(cv){
        ets_mae <- abs(as.vector(ets_base$mean) - as.vector(y_test))
        ets_base <- forecast::forecast(forecast::ets(y),
                                       h = frequency * k)
      } else {
        ets_mae <- tail(as.vector(abs(residuals(ets_base))), length(y_test))
      }
    }
  } else {
    ets_base <- try(forecast::forecast(forecast::tbats(y_train),
                                       h = length(y_test)),
                    silent = TRUE)
    if(inherits(ets_base, 'try-error')){
      use_ets <- FALSE
      ets_mae <- rep(NA, length(y_test))
    } else {
      use_ets <- TRUE
      if(cv){
        ets_mae <- abs(as.vector(ets_base$mean) - as.vector(y_test))
        ets_base <- forecast::forecast(forecast::tbats(y),
                                       h = frequency * k)
      } else {
        ets_mae <- tail(as.vector(abs(residuals(ets_base))), length(y_test))
      }
    }
  }

if(!bottom_series){
# Try an auto.arima model
arima_base <- try(forecast::forecast(forecast::auto.arima(y_train,
                                                          lambda = lambda),
                                     h = length(y_test)),
                                     silent = TRUE)

if(inherits(arima_base, 'try-error')){
  use_arima <- FALSE
  arima_mae <- rep(NA, length(y_test))
} else {
  use_arima <- TRUE
  if(cv){
    arima_mae <- abs(as.vector(arima_base$mean) - as.vector(y_test))

    arima_base <- forecast::forecast(forecast::auto.arima(y_train,
                                                          lambda = lambda),
                                     h = frequency * k)
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
  maes <- cbind(stlf_mae,
                ets_mae,
                arima_mae,
                arimaf_mae,
                naive_mae,
                rwf_mae,
                snaive_mae)
  colnames(maes) <- c('stlf', 'ets', 'arima', 'arimaf',
                      'naive', 'rwf', 'snaive')
} else {
  maes <- cbind(stlf_mae,
                ets_mae,
                naive_mae,
                rwf_mae,
                snaive_mae)
  colnames(maes) <- c('stlf', 'ets', 'naive', 'rwf', 'snaive')
}


if(any(which(apply(maes, 2, FUN = function(x) length(which(is.na(x)))) == NROW(maes)))){
  maes <- maes[,-which(apply(maes, 2, FUN = function(x)
    length(which(is.na(x)))) == NROW(maes))]
}

# Scale against the snaive model to convert the MAE to MASE
# The scaling should depend on seasonality and should use naive when there is no seasonality.
# But experimentation suggests the snaive will be
# equivalent to the naive when seasonality is not present
maes <- maes + 1 / (snaive_mae + 1)

# Define the function to minimise
opt_mae <- function(weights, maes) {
  mean.default(unlist(lapply(seq_len(NROW(maes)), function(h){
    weighted.mean(maes[h,], weights)
  }), use.names = FALSE))
}

# Use optimisation to minimise the function
weights <- rbeta(NCOL(maes), 1, 1)
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
  ens_mean <- unlist(lapply(seq_len(NROW(means)), function(x){
    weighted.mean(means[x,], ens_weights)}), use.names = F)
  rm(means)

  upper_95s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$upper[,2])}))
  ens_upper95 <- unlist(lapply(seq_len(NROW(upper_95s)), function(x){
    weighted.mean(upper_95s[x,], ens_weights)}), use.names = F)
  rm(upper_95s)

  upper_80s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$upper[,1])}))
  ens_upper80 <- unlist(lapply(seq_len(NROW(upper_80s)), function(x){
    weighted.mean(upper_80s[x,], ens_weights)}), use.names = F)
  rm(upper_80s)

  lower_95s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$lower[,2])}))
  ens_lower95 <- unlist(lapply(seq_len(NROW(lower_95s)), function(x){
    weighted.mean(lower_95s[x,], ens_weights)}), use.names = F)
  rm(lower_95s)

  lower_80s <- do.call(cbind, lapply(seq_along(fcs), function(x){
    as.vector(fcs[[x]]$lower[,1])}))
  ens_lower80 <- unlist(lapply(seq_len(NROW(lower_80s)), function(x){
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
  ens_fitted <- unlist(lapply(seq_len(NROW(fitted)), function(x){
    weighted.mean(fitted[x,], ens_weights)}), use.names = F)
  residuals <- ens_fitted - y
  residuals[is.na(residuals)] <- 0

  return(list(ensemble_forecast, residuals))
}

# Gather the base forecasts and calculate the ensemble
if(!bottom_series){
  fcs <- list(stlf_base,
              ets_base,
              arima_base,
              arimaf_base,
              naive_base,
              rwf_base,
              snaive_base)
  names(fcs) <- c('stlf', 'ets', 'arima', 'arimaf',
                  'naive', 'rwf', 'snaive')
} else {
  fcs <- list(stlf_base,
              ets_base,
              naive_base,
              rwf_base,
              snaive_base)
  names(fcs) <- c('stlf', 'ets', 'naive', 'rwf', 'snaive')
}

fcs <- fcs[names(fcs) %in% names(ens_weights)]
ensemble <- weight_fcs(fcs = fcs, ens_weights = ens_weights)
rm(fcs, ens_weights, maes,
   ets_base,
   stlf_base,
   naive_base,
   rwf_base,
   snaive_base)

# Return the ensemble forecast and residuals
if(discrete){

  # Transform the forecast and residuals to the estimated discrete distribution
  trans_fc <- transform_fc_preds(.subset2(ensemble, 1),
                                 copula_y$y_discrete, dist_mappings, dist_params)
  if(length(dist_params) == 2){
    trans_resids <- ts(back_trans(as.vector(.subset2(ensemble, 2)),
                                  dist_mappings, dist_params) - dist_params[2],
                       start = start(.subset2(ensemble, 2)),
                       frequency = frequency(.subset2(ensemble, 2)))
  } else {
    trans_resids <- ts(back_trans(as.vector(.subset2(ensemble, 2)),
                                  dist_mappings, dist_params) - dist_params,
                       start = start(.subset2(ensemble, 2)),
                       frequency = frequency(.subset2(ensemble, 2)))
  }
  output <- list(forecast = trans_fc, residuals = trans_resids)

} else {
  output <- list(forecast = .subset2(ensemble, 1), residuals = .subset2(ensemble, 2))
}
return(output)
}
