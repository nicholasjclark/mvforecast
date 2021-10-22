#'Fit an stlf model to a logged version of a discrete time series
#'
#'This function fits an stlf model with AR errors to a discrete time series that has been interpolated
#'and log(x + 1/6) transformed
#'
#'
#'@param y \code{xts vector}. The discrete outcome series to be modelled. \code{NAs} are allowed
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param horizon \code{integer}. The horizon to forecast. Defaults to \code{frequency}
#'
#'@return A \code{vector} containing fitted values from the model of \code{length(length(y)+horizon)}
#'@export
#'
fit_stlf_discreet = function(y, frequency, horizon){

  # Check variables
  if (!xts::is.xts(y)) {
    stop("y must be an xts object")
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

  # Impute any NAs
  y_impute <- round(as.vector(forecast::na.interp(zoo::na.approx(zoo::na.approx(ts(
    y, frequency = frequency), na.rm = F), na.rm = F))), 0)

  y_impute <- ts(y_impute, frequency = frequency,
                 start = start(xts.to.ts(y, freq = frequency)))
  y_impute[y_impute < 0] <- 0

  # Log(y + 1/6) for ARIMA modelling
  log_y <- ts(log(y_impute + 1/6), start = start(y_impute),
              frequency = frequency(y_impute))

  # Fit auto.arima()
  fc_mod <- forecast::stlm(log_y, modelfunction = ar)
  fc_fc <- forecast::forecast(fc_mod, horizon)

  # Extract fitted and forecasted values
  as.vector(forecast::na.interp(zoo::na.approx(c(fc_mod$fitted, fc_fc$mean), na.rm = F)))
}
