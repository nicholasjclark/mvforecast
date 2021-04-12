#'Test script to ensure mvforecast is working properly
#'
#'This function fits a vets model and equivalent thief_vets model on the example ixodes_vets_dat dataset
#'and returns a plot and a CRPS calculation
#'
#'@param plot \code{logical}. Default is \code{TRUE}
#'@return A plot from \code{\link{plot_vets_preds}} and a CRPS calculation from \code{\link{calc_crps}}
#'@export
test_mvforecast = function(plot = TRUE){
  cat('Loading the ixodes_vets_dat dataset\n\n')
  data("ixodes_vets_dat")

  cat('Fitting a vets model with common seasonal smoothing and grouped level smoothing\n\n')
  mod1 <- tsvets:::simulate.tsvets.estimate(tsvets:::estimate.tsvets.spec(tsvets:::vets_modelspec(ixodes_vets_dat$y_train,
                                           level = "grouped",
                                           slope = "none",
                                           damped = "none",
                                           seasonal = "common",
                                           lambda = 0.75,
                                           dependence = "equicorrelation",
                                           frequency = 52,
                                           cores = parallel::detectCores() - 1,
                                           group = ixodes_vets_dat$groups),
                            solver = "solnp",
                            control = list(trace = 0)),
                   nsim = 1000,
                   h = ixodes_vets_dat$h)

  cat('Fitting an equivalent model that then uses hierarchical reconciliation\n\n')
  mod2 <- thief_vets(y = ixodes_vets_dat$y_train,
                     multi_freq = 30,
                     level = "grouped",
                     slope = "none",
                     damped = "none",
                     seasonal = "common",
                     lambda = 0.75,
                     dependence = "equicorrelation",
                     frequency = 52,
                     cores = parallel::detectCores() - 1,
                     group = ixodes_vets_dat$groups,
                     save_plots = F)

  if(plot){
    cat('\n\nPlotting thief simulation forecast (lines) and true values (ytest points) for NEON plot_ID 4\n\n')
    plot_mv_preds(simulation = mod2[[4]])
    points(as.vector(ixodes_vets_dat$y_test[,4]))
  }

  cat('\n\nCalculating CRPS against ixodes_vets_dat$ytest for both models (lower is better)\n\n')
  vets_crps <- calc_crps(mod1, y_test = ixodes_vets_dat$y_test)
  thief_vets_crps <- calc_crps(simulation = mod2, y_test = ixodes_vets_dat$y_test)

  return(data.frame(rbind(vets_crps, thief_vets_crps)))
}

