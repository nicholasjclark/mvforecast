#'Probabilistic reconciliation of temporally reconciled grouped forecasts
#'
#'This function fits a specified thief model and then uses probabilistic reconciliation based on the supplied
#'grouping structure to reconcile base forecast distributions
#'
#'@importFrom parallel detectCores parLapply makePSOCKcluster setDefaultCluster clusterExport clusterEvalQ stopCluster
#'@importFrom stats ts end start frequency
#'
#'@param y \code{xts matrix}. The outcome series to be modelled. \code{NAs} are currently not supported
#'@param k \code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}. Default
#'is \code{1}, meaning that a final forecast of \code{frequency} horizons will be returned
#'@param lambda \code{numeric}. The Box Cox power transformation parameter for aggregate series. Must be
#'between \code{-1} and \code{2} inclusive. If \code{y_series} contains zeros, \code{lambda} will be set to
#'\code{max(c(0.7, lambda))} to ensure stability of forecasts
#'@param model One of code{c(thief_ensemble, thief_rfsrc, thief_vets)}
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param horizon \code{integer}. The horizon to forecast. Defaults to \code{frequency}
#'@param groups Group matrix indicates the group structure, with one column for each series when
#'completely disaggregated, and one row or each grouping of the time series. It allows either a numerical matrix or a
#'matrix consisting of strings that can be used for labelling. See \code{\link[hts]{gts}} for more details
#'@param cores \code{integer}. The number of cores to use
#'@param keep_total \code{logical}. If \code{TRUE}, forecasts for the the top-level summed series (the total)
#'will also be returned as the last element in the returned list of forecast distributions
#'@param max_agg (optional) \code{integer} specifying the maximum number of temporal aggregation levels
#'to use when reconciling, via the structural scaling method. Useful if higher levels of aggregation
#'are unlikely to have 'seen' recent changes in series dynamics and will likely then result in poor
#'forecasts as a result. Default is \code{NULL}, meaning that all levels of aggregation are used
#'@param ... Other arguments to be passed on to the specified thief model
#'@return A \code{list} containing the reconciled forecast distributions for each series in \code{y}. Each element in
#'the \code{list} is a \code{horizon x 1000 matrix} of forecast predictions
#'
#'@seealso \code{\link{thief_ensemble}}, \code{\link[hts]{gts}},
#'\code{\link[thief]{reconcilethief}}, \code{\link[ProbReco]{scoreopt}}
#'
#'@details Series in \code{y} are aggregated at all possible levels up to annual using \code{\link[thief]{tsaggregates}}.
#'The specified model is used on the unaggregated series, with forecasts reconciled
#'using \code{\link[thief]{reconcilethief}}. Finally, the grouping structure is used in \code{\link[hts]{gts}} to create
#'the group summing structure. Forecasts of the summed series, including the final total, are generated using
#'\code{\link{thief_ensemble}} and the final base forecast distributions are reconciled using weights calculated through
#'energy score optimisation using \code{\link[ProbReco]{scoreopt}}
#'
#'@references Athanasopoulos, G., Hyndman, R.,  Kourentzes, N., and Petropoulos, F. Forecasting with temporal hierarchies.
#'(2017) European Journal of Operational Research 262(1) 60–74 \cr\cr
#'Panagiotelis, A., Gamakumara, P., Athanasopoulos, G., and Hyndman, R.
#'Probabilistic Forecast Reconciliation: Propoerties, Evaluation and Score Optimisation. (2020)
#'Monash EBS Working Paper 26/20 \cr\cr
#'Panagiotelis, A., Athanasopoulos, G., Gamakumara, P., and Hyndman, R. (2021). Forecast reconciliation:
#'A geometric view with new insights on bias correction. International Journal of Forecasting, 37(1), 343-359.
#'
#'@examples
#'\donttest{
#'library(mvforecast)
#'data("ixodes_vets_dat")
#'
#'#Fit a probabilistic thief_ensemble model
#'groups <- ixodes_vets_dat$groups
#'mod1 <- thief_probreco(y = ixodes_vets_dat$y_train,
#'model = 'thief_ensemble',
#'frequency = 52,
#'prob_train_horizon = 12,
#'keep_total = FALSE,
#'groups = groups)
#'
#'#Calculate the out-of-sample CRPS
#'calc_crps(mod1, y_test = ixodes_vets_dat$y_test)
#'
#'Plot simulation results for one of the plots in the NEON dataset
#'plot_mvforecast(simulation = mod1[[4]])
#'points(as.vector(ixodes_vets_dat$y_test[,4]))}
#'
#'@export
#'
thief_probreco = function(y,
                          k = 1,
                          frequency = 52,
                          lambda = NULL,
                          model = 'thief_ensemble',
                          prob_train_horizon = NULL,
                          horizon = NULL,
                          max_agg = NULL,
                          groups,
                          keep_total = TRUE,
                          cores = parallel::detectCores() - 1, ...){

# Check variables
  if (!xts::is.xts(y)) {
    stop("y must be an xts object")
  }

  model <- match.arg(arg = model[1],
                     choices = c("thief_ensemble", "thief_rfsrc", "thief_vets"))

# Function to convert xts to ts object
  xts.to.ts <- function(x, freq = 52) {
    start_time <- floor((lubridate::yday(start(x)) / 365) * freq)
    ts(as.numeric(x),
       start = c(lubridate::year(start(x)),
                 start_time), freq = freq)
  }

# Split y so that last prob_train_horizon is used for training probabilistic weights
original_y <- y
y_test <- y[(nrow(y) - prob_train_horizon + 1):nrow(y), ]
y <- y[1:(nrow(y) - prob_train_horizon), ]

# Fit the multivariate reconciled model for the base training series
if(model == 'thief_rfsrc'){
  cat('Fitting thief_rfsrc for the base series\n\n')
  y_forecast_probs <- thief_rfsrc(y = y,
                                  frequency = frequency,
                                  k = ceiling(prob_train_horizon / frequency),
                                  horizon = prob_train_horizon,
                                  cores = cores,
                                  max_agg = max_agg,
                                  ...)
}

if(model == 'thief_vets'){
  cat('Fitting thief_vets for the base series\n\n')
  y_forecast_probs <- try(thief_vets(y = y,
                                     lambda = lambda,
                                     k = ceiling(prob_train_horizon / frequency),
                                     frequency = frequency,
                                     horizon = prob_train_horizon,
                                     cores = cores,
                                     max_agg = max_agg,
                                     group = groups, ...), silent = T)

  # Errors will occur if there aren't enough seasonal periods in the data for ETS
  # Use thief_rfsrc if this happens, as this is still a multivariate base model
  if(inherits(y_forecast_probs, 'try-error')){
    cat('**thief_vets errored due to too few seasonal periods in training data.\nUsing thief_rfsrc instead**\n\n')
    model <- 'thief_rfsrc'
    y_forecast_probs <- thief_rfsrc(y = y,
                                    frequency = frequency,
                                    k = ceiling(prob_train_horizon / frequency),
                                    horizon = prob_train_horizon,
                                    cores = cores,
                                    max_agg = max_agg,
                                    ...)
  }
}

if(model == 'thief_ensemble'){
  cat('Fitting thief_ensemble for the base series\n\n')
  y_forecast_probs <- thief_ensemble(y = y,
                                     frequency = frequency,
                                     k = ceiling(prob_train_horizon / frequency),
                                     cores = cores,
                                     max_agg = max_agg,
                                     horizon = prob_train_horizon)
}

# Calculate aggregated group series, including the final total
y_ts <- do.call(cbind,lapply(seq_len(ncol(y)), function(x){
  xts.to.ts(y[,x], frequency)
}))
colnames(y_ts) <- colnames(y)
train <- hts::gts(y_ts, groups = t(matrix(groups)), gnames = 'species')
aggregated_series <- hts::aggts(train)[,1:(ncol(hts::aggts(train)) - ncol(y_ts))]
S <- hts::smatrix(train)

# Fit thief_ensemble models to the aggregated series
cat('\n\nFitting thief_ensemble for the aggregated group series\n\n')
all_agg_series <- zoo::as.zoo(aggregated_series)
all_agg_series <- xts::xts(all_agg_series, lubridate::date_decimal(zoo::index(all_agg_series)))
group_reconciled <- thief_ensemble(y = all_agg_series,
                                   frequency = frequency,
                                   cores = cores,
                                   k = ceiling(prob_train_horizon / frequency),
                                   max_agg = max_agg,
                                   horizon = prob_train_horizon)

# Get realisations from the out of sample test set for all series
y_ts_test <- do.call(cbind,lapply(seq_len(ncol(y_test)), function(x){
  xts.to.ts(y_test[,x], frequency)
}))
colnames(y_ts_test) <- colnames(y_test)
test <- hts::gts(y_ts_test, groups = t(matrix(groups)), gnames = 'species')
aggregated_series_test <- hts::aggts(test)[,1:(ncol(hts::aggts(test)) - ncol(y_ts))]
data_matrix <- rbind(t(aggregated_series_test[1:nrow(y_test),]),
                     as.matrix(t(y_test)))
reco_data <- lapply(seq_len(NCOL(data_matrix)), function(y){
  matrix(data_matrix[,y])
})

# Gather probabilistic distributions for base and grouped series
prob_distributions <- lapply(seq_len(nrow(y_test)), function(h){
  rbind(do.call(rbind, lapply(seq_along(group_reconciled), function(x){
    t(group_reconciled[[x]][h, ])
  })),
  do.call(rbind, lapply(seq_along(y_forecast_probs), function(series){
    t(y_forecast_probs[[series]][h, ])
  })))
})

# Function needed to simulate from each forecast distribution
make_genfunc <- function(input){
  f <- function(){
    out <- matrix(0, n_series, 1000)
    for (j in 1:n_series){
      out[j,] <- sample(input[j,], 1000, replace = T)
    }

    return(out)
  }
  return(f)
}

# Map the function across each element in prob_distributions
reco_prob <- purrr::map(prob_distributions, make_genfunc)

# Optimise reconciliation weigths to minimise prediction interval error
cat('\n\nOptimising out-of-sample prediction interval reconciliation using ProbReco\n\n')
n_series <- ncol(hts::aggts(train))
opt <- ProbReco::scoreopt(reco_data,
                          reco_prob,
                          S,
                          score = list(score = 'energy', alpha = 1),
                          control = list(maxIter = 500, tol = 1E-08),
                          matches = TRUE)

# Generate distributions for the out-of-sample testing period
cat('Optimised weights calculated. Re-calibrating models on the full supplied y data\n\n')
# Set forecast horizon if missing
if(missing(horizon)){
  horizon <- frequency
}

if(model == 'thief_rfsrc'){
  cat('Fitting thief_rfsrc for the full base series\n\n')
  y_forecast_probs <- thief_rfsrc(y = original_y,
                                  frequency = frequency,
                                  k = ceiling(horizon / frequency),
                                  horizon = horizon,
                                  cores = cores,
                                  max_agg = max_agg,
                                  ...)
}

if(model == 'thief_vets'){
  cat('Fitting thief_vets for the full base series\n\n')
  y_forecast_probs <- try(thief_vets(y = original_y,
                                     lambda = lambda,
                                     k = ceiling(horizon / frequency),
                                     frequency = frequency,
                                     horizon = horizon,
                                     cores = cores,
                                     max_agg = max_agg,
                                     group = groups, ...), silent = T)

  # Errors will occur if there aren't enough seasonal periods in the data for ETS
  # Use thief_rfsrc if this happens, as this is still a multivariate base model
  if(inherits(y_forecast_probs, 'try-error')){
    cat('**thief_vets errored due to too few seasonal periods in training data.\nUsing thief_rfsrc instead**\n\n')
    model <- 'thief_rfsrc'
    y_forecast_probs <- thief_rfsrc(y = original_y,
                                    cores = cores,
                                    frequency = frequency,
                                    k = ceiling(horizon / frequency),
                                    max_agg = max_agg,
                                    horizon = horizon, ...)
  }
}

if(model == 'thief_ensemble'){
  cat('Fitting thief_ensemble for the full base series\n\n')
  y_forecast_probs <- thief_ensemble(y = original_y,
                                     cores = cores,
                                     frequency = frequency,
                                     k = ceiling(horizon / frequency),
                                     max_agg = max_agg,
                                     horizon = horizon)
}

# Generate aggregates for the full supplied y data
y_ts <- do.call(cbind,lapply(seq_len(ncol(original_y)), function(x){
  xts.to.ts(original_y[,x], frequency)
}))
colnames(y_ts) <- colnames(original_y)
orig <- hts::gts(y_ts, groups = t(matrix(groups)), gnames = 'species')

aggregated_series <- hts::aggts(orig)[,1:(ncol(hts::aggts(orig)) - ncol(y_ts))]

# Fit thief_ensemble models to the full aggregated series
cat('\n\nFitting thief_ensemble for the full aggregated group series\n\n')
all_agg_series <- zoo::as.zoo(aggregated_series)
all_agg_series <- xts::xts(all_agg_series, lubridate::date_decimal(zoo::index(all_agg_series)))
group_reconciled <- thief_ensemble(y = all_agg_series,
                                   frequency = frequency,
                                   k = ceiling(horizon / frequency),
                                   horizon = horizon,
                                   max_agg = max_agg,
                                   cores = cores)

# Extract future forecast distributions for all series
new_prob_distributions <- lapply(seq_len(nrow(y_forecast_probs[[1]])), function(h){
  rbind(do.call(rbind, lapply(seq_along(group_reconciled), function(x){
    t(group_reconciled[[x]][h, ])
  })),
  do.call(rbind, lapply(seq_along(y_forecast_probs), function(series){
    t(y_forecast_probs[[series]][h, ])
  })))
})

# Reconcile forecasts using weights from ProbReco
optimised_base <- lapply(seq_len(nrow(y_forecast_probs[[1]])), function(h){
  base_opt <- S %*% (opt$d + opt$G %*% new_prob_distributions[[h]])
  # Drop intermediate group totals, but keep the final total
  if(keep_total){
    order <- c((ncol(aggregated_series) + 1):nrow(base_opt),1)
    base_opt <- base_opt[order,]
  } else {
    base_opt <- base_opt[-c(1: ncol(aggregated_series_test)),]
  }

  if(!any(y < 0)){
    base_opt[base_opt < 0] <- 0
  }
  base_opt
})

# Put reconciled forecast distributions into correct format
reconciled_base <- lapply(seq_len(nrow(optimised_base[[1]])), function(x){
  all <- do.call(rbind, lapply(seq_along(optimised_base), function(y){
    optimised_base[[y]][x,]
  }))
  all
})

# Return the reconciled forecast distributions for each series in y as a list
if(keep_total){
  names(reconciled_base) <- c(colnames(y), 'total')
} else {
  names(reconciled_base) <- colnames(y)
  }

return(reconciled_base)
}
