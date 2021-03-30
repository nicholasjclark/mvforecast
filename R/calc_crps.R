#'calculate continuous rank probability score from a VETS model forecast distribution
#'
#'
#'@importFrom parallel detectCores
#'@importFrom stats ts end start frequency
#'
#'@param simulation \code{matrix}. The forecast object, which can either be a simulation from
#'\code{\link[tsvets]{predict.tsvets.estimate}} or a list of reconciled distributions from \code{\link{thief_vets}}
#'@param y_test \code{xts matrix}. The out-of-sample test horizon of the outcome series
#'
#'@details The sum of the CRPS is calculated against the test observations in \code{y_test} and returned. Calculations
#'use the \code{\link[scoringRules]{crps_sample}} function
#'
#'@return a \code{summary} of the out-of-sample CRPS
#'
#'@export
calc_crps = function(simulation, y_test){
  crps_sums <- vector()

  # Check class of simulation
  if(class(simulation) == "tsvets.simulate"){
    raw_distribution <- FALSE
  } else {
    raw_distribution <- TRUE
  }

  # Calculate sum of the CRPS across the test horizon specified by y_test
  if(raw_distribution){
    for(i in 1:ncol(y_test)){
      matched_preds <- simulation[[i]][1:nrow(y_test),]
      crps_sums[i] <- sum(scoringRules::crps_sample(y = as.vector(y_test[,i]),
                                                    dat = matched_preds))
    }
  } else {
    # If this is a vets prediction object, get the predictions from the correct list slot
    for(i in 1:ncol(y_test)){
      simulation$simulation_table$Simulated[[i]]$distribution[
        simulation$simulation_table$Simulated[[i]]$distribution < 0] <-0
      crps_sums[i] <- sum(scoringRules::crps_sample(y = as.vector(y_test[,i]),
                                                    dat = t(simulation$simulation_table$Simulated[[i]]$distribution)))
    }

  }
  summary(crps_sums)
}
