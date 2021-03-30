#'Test script to ensure mvforecast is working properly
#'
#'This function fits a thief_vets model on the example ixodes_vets_dat dataset and returns a plot and a CRPS calculation
#'
#'@param plot \code{logical}. Default is \code{TRUE}
#'@return A plot from \code{\link{plot_vets_preds}} and a CRPS calculation from \code{\link{calc_crps}}
#'@export
test_mvforecast = function(plot = TRUE){
  cat('Loading the ixodes_vets_dat dataset\n\n')
  data("ixodes_vets_dat")
  mod <- thief_vets(y = ixodes_vets_dat$y_train,
                     multi_freq = 30,
                     level = "grouped",
                     slope = "none",
                     damped = "none",
                     seasonal = "common",
                     lambda = 1,
                     dependence = "equicorrelation",
                     frequency = 52,
                     cores = parallel::detectCores() - 1,
                     group = ixodes_vets_dat$groups,
                     save_plots = F)

  if(plot){
    cat('\n\nPlotting simulation forecast (lines) and true values (ytest points) for NEON plot_ID 4\n\n')
    plot_vets_preds(simulation = mod[[4]])
    points(as.vector(ixodes_vets_dat$y_test[,4]))
  }

  cat('\n\nCalculating CRPS against ixodes_vets_dat$ytest\n\n')
  return(calc_crps(simulation = mod, y_test = ixodes_vets_dat$y_test))
}

