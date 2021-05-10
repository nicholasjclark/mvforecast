#'Plot simulation predictions from a reconciled thief model
#'
#'
#'
#'@param simulation \code{matrix}. The forecast object, which should be a horizon x simulation \code{matrix}
#'similar to those returned by \code{\link{thief_vets}}
#'@param main \code{character}. Optional label for the plot
#'@param ylab \code{character}. Optional label for the y axis
#'@param ylim Limits for the y-axis. Default is \code{c(min(simulation, na.rm = T), max(simulation, na.rm = T))}
#'
#'@details The forecast distribution is used to calculate highest posterior density credible estimates for the following
#'intervals: \code{c(0.95, 0.9, 0.85, 0.8, 0.75, 0.7)}
#'
#'@return A base plot of the forecast medeian (red dashed line) and credible intervals as a gradient polygon
#'
#'@export
plot_mvforecast = function(simulation, main = '', ylab = 'Y',
                         ylim = c(min(simulation, na.rm = T), max(simulation, na.rm = T))){

  # Calculate prediction intervals
  forecast <- do.call(rbind, lapply(seq_len(nrow(simulation)), function(x){
    pred_vals <- simulation[x, ]
    pred_vals <- pred_vals[!is.na(pred_vals)]
    ninetyfives <- suppressWarnings(hpd(pred_vals, 0.95))
    nineties <- suppressWarnings(hpd(pred_vals, 0.9))
    eightyfives <- suppressWarnings(hpd(pred_vals, 0.85))
    eighties <- suppressWarnings(hpd(pred_vals, 0.8))
    seventyfives <- suppressWarnings(hpd(pred_vals, 0.75))
    seventies <- suppressWarnings(hpd(pred_vals, 0.70))
    quantiles <- c(ninetyfives[1], nineties[1], eightyfives[1], eighties[1], seventyfives[1], seventies[1],
                   eighties[2],
                   seventies[3], seventyfives[3], eighties[3], eightyfives[3], nineties[3], ninetyfives[3])
    quantiles
  }))

  # Create an empty plot
  if(nrow(forecast) > 20){
    x_break_length <- 5
  } else if(nrow(forecast) > 10){
    x_break_length <- 4
  } else if(nrow(forecast) > 5){
    x_break_length <- 2
  } else {
    x_break_length <- 1
  }

  if(max(forecast, na.rm = T) > 10){
    y_break_length <- 6
  } else {
    y_break_length <- 4
  }

  plot(as.vector(simulation[,1]), xlab = 'Horizon', main = main,
       ylab = ylab,
       ylim = c(min(simulation, na.rm = T), max(forecast, na.rm = T)), type = 'n',
       xaxt='n',yaxt='n', axes = FALSE)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(0.98, 0.98, 0.98, 1),
       border = NA)
  axis(side = 2, at = round(seq(0, max(forecast, na.rm = T), length.out = y_break_length), 1))
  axis(side = 1, at = seq(0, nrow(forecast) + x_break_length, by = x_break_length))


  # Add mean and 95% HPD intervals as lines
  polygon(c(seq(1:nrow(forecast)), rev(seq(1:nrow(forecast)))),
          c(forecast[,1] , rev(forecast[,13])), col = rgb(1, 0, 0, 0.075),
          border = NA)
  polygon(c(seq(1:nrow(forecast)), rev(seq(1:nrow(forecast)))),
          c(forecast[,2] , rev(forecast[,12])), col = rgb(1, 0, 0, 0.1),
          border = NA)
  polygon(c(seq(1:nrow(forecast)), rev(seq(1:nrow(forecast)))),
          c(forecast[,3] , rev(forecast[,11])), col = rgb(1, 0, 0, 0.125),
          border = NA)
  polygon(c(seq(1:nrow(forecast)), rev(seq(1:nrow(forecast)))),
          c(forecast[,4] , rev(forecast[,10])), col = rgb(1, 0, 0, 0.15),
          border = NA)
  polygon(c(seq(1:nrow(forecast)), rev(seq(1:nrow(forecast)))),
          c(forecast[,5] , rev(forecast[,9])), col = rgb(1, 0, 0, 0.175),
          border = NA)
  polygon(c(seq(1:nrow(forecast)), rev(seq(1:nrow(forecast)))),
          c(forecast[,6] , rev(forecast[,8])), col = rgb(1, 0, 0, 0.2),
          border = NA)
  lines(forecast[,7], col = rgb(1, 0, 0, 0.75), lwd = 2, lty = 'dashed')
  for(i in seq(0, nrow(forecast), by = x_break_length/2)){
    abline(v=i, col = rgb(0.98, 0.98, 0.98, 0.45), lwd = 0.05)
  }
  for(i in seq(0, max(forecast, na.rm = T), length.out = y_break_length*2)){
    abline(h=i, col = rgb(0.98, 0.98, 0.98, 0.45), lwd = 0.05)
  }

}
