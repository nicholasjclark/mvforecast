#'Plot simulation predictions from a reconciled VETS model
#'
#'
#'
#'@param simulation \code{matrix}. The forecast object, which should be a horizon x simulation \code{matrix}
#'similar to those returned by \code{\link{thief_vets}}
#'@param main \code{character}. Optional label for the plot
#'@param ylab \code{character}. Optional label for the y axis
#'
#'@details The sum of the CRPS is calculated against the test observations in \code{y_test} and returned. Calculations
#'use the \code{\link[scoringRules]{crps_sample}} function
#'
#'@return A base plot of the forecast mean (black line) and 95 percent prediction intervals (red dashed lines)
#'
#'@export
plot_mv_preds = function(simulation, main = '', ylab = 'Y'){

  # Create an empty plot
  plot(as.vector(simulation[,1]), xlab = 'Horizon', main = main,
       ylab = ylab,
       ylim = c(min(simulation, na.rm = T), max(simulation, na.rm = T)), type = 'n')

  # Calculate prediction intervals
  forecast <- do.call(rbind, lapply(seq_len(nrow(simulation)), function(x){
    pred_vals <- simulation[x, ]
    pred_vals <- pred_vals[!is.na(pred_vals)]
    ninetyfives <- suppressWarnings(hpd(pred_vals, 0.95))
    eighties <- suppressWarnings(hpd(pred_vals, 0.8))
    quantiles <- c(ninetyfives[1], eighties[1], mean(pred_vals, na.rm = T), eighties[3], ninetyfives[3])
    quantiles
  }))

  # Add mean and 95% HPD intervals as lines
  lines(forecast[,3], col = 'black', lwd = 2)
  lines(forecast[,1], col = 'red', lty = 2)
  lines(forecast[,5], col = 'red', lty = 2)
}
