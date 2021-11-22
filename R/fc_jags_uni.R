#'Fit a Bayesian forecast model to a univariate discrete time series
#'
#'This function takes fitted values from a log-transformed forecast model and uses those as
#'a predictor in a Negative Binomial regression to predict the discrete outcome time series
#'
#'
#'@param y \code{xts vector}. The discrete outcome series to be modelled. \code{NAs} are allowed
#'@param model \code{character} either 'TBATS', 'arima' or 'stlf'
#'@param frequency \code{integer}. The seasonal frequency in \code{y}
#'@param horizon \code{integer}. The horizon to forecast. Defaults to \code{frequency}
#'@param use_nb \code{logical} If \code{TRUR}, use a Negative Binomial likelihood with estimated
#'overdispersion parameter r;
#'@param n.chains \code{integer} the number of parallel chains for the model
#'@param n.adapt \code{integer} the number of iterations for adaptation. See adapt for details.
#'If \code{n.adapt} = 0 then no adaptation takes place
#'@param n.iter \code{integer} the number of iterations of the Markov chain to run
#'@param auto_update \code{logical} If \code{TRUE}, the model is run for up to \code{3} additional sets of
#'\code{n.iter} iterations, or until the lower 15th percentile of effective sample sizes reaches \code{200}
#'
#'@return A \code{list} containing the posterior forecast from the forecast model and the original fitted
#'\code{\link[rjags]{jags.model}}
#'
#'@details The discrete series is first log-transformed and interpolated to fit a forecast model. Fitted values
#'from the forecast model are then used as a predictor in a Negative Binomial Bayesian regression to model
#'the raw discrete series
#'
#'@export
#'
fc_jags_uni = function(y,
                       model = 'arima',
                       horizon,
                       frequency,
                       use_nb = TRUE,
                       n.chains = 2,
                        n.adapt = 1000,
                        n.burnin = 10000,
                          n.iter = 10000,
                          thin = 10,
                          auto_update = TRUE){
if(use_nb){
  use_nb = 1
} else {
  use_nb = 0
}

if(!tolower(model) %in% c('tbats', 'arima', 'stlf')){
  stop("model must be either 'TBATS', 'arima' or 'stlf'")
}

  if(tolower(model) == 'tbats'){
    fc_fit <- fit_tbats_discreet(y,
                                 frequency = frequency,
                                 horizon = horizon)
  }

  if(tolower(model) == 'arima'){
    fc_fit <- fit_arima_discreet(y,
                                 frequency = frequency,
                                 horizon = horizon)
  }

  if(tolower(model) == 'stlf'){
    fc_fit <- fit_stlf_discreet(y,
                                 frequency = frequency,
                                 horizon = horizon)
  }

y_orig <- y
y <- c(as.vector(y), rep(NA, horizon))

model_file <- " model {

## Mean expectations
 for (i in 1:n) {
   mu[i] <- exp(alpha + fc_fit[i])
 }
# alpha should be negative as we added 1/6 prior to the univariate model
alpha ~ dunif(-2, 0)

## Negative binomial likelihood
for (i in 1:n) {
 y[i] ~ dnegbin(rate[i], r);
 rate[i] <- ifelse((r / (r + mu[i])) < min_eps, min_eps,
                  (r / (r + mu[i])))
}

## Set overdispersion parameter to 10000 if not using nb;
## this approximates the Poisson distribution
 r1 ~ dgamma(0.01, 0.01)
 r2 <- 10000
 r_indicator <- ifelse(use_nb > 0, 1, 0)
 r <- (r1 * r_indicator) + (r2 * (1 - r_indicator))

 ## Posterior predictions
 for (i in 1:n) {
  ypred[i] ~ dnegbin(rate[i], r)
 }

}"

jags_data <- list(
  y = y,
  fc_fit = fc_fit,
  n = length(y),
  min_eps = .Machine$double.eps,
  use_nb = use_nb
)

inits <- function(){
  list(alpha = -0.5)
}

# Initiate adaptation of the model
load.module("glm")
jags_mod <- jags.model(textConnection(model_file),
                      data = jags_data,
                      inits = inits,
                      n.chains = n.chains,
                      n.adapt = 0)


# Update the model for the burnin period
adapt(jags_mod, n.burnin, end.adaptation = TRUE)

# Sample from the posterior
out_jags_mod <- rjags::coda.samples(jags_mod,
                              variable.names = c('alpha',
                                                 'ypred',
                                                 'r'),
                              n.iter = n.iter,
                              thin = thin)

if(auto_update){
  # Update until reasonable convergence in the form of Rhat and ESS
    update_params <- c('alpha', 'r')

  mod_summary <- MCMCvis::MCMCsummary(out_jags_mod, update_params)
  for(i in 1:3){
    if(quantile(mod_summary$Rhat, 0.85, na.rm = T) <= 1.1 &
       quantile(mod_summary$n.eff, 0.15) >= 200){
      break;
    }
    cat('Convergence not reached. Extending burnin...\n')
    update(jags_mod, n.burnin)
    out_jags_mod <- rjags::coda.samples(jags_mod,
                                       variable.names = c('alpha',
                                                          'ypred',
                                                          'r'),
                                       n.iter = n.iter,
                                       thin = thin)
    mod_summary <- MCMCvis::MCMCsummary(out_jags_mod, update_params)
  }
}

# Return forecast and fitted jags model
 ypreds <- t(MCMCvis::MCMCchains(out_jags_mod, 'ypred'))
 ypreds <- ypreds[(length(y_orig)+1):NROW(ypreds), ]
 list(forecast = ypreds, jags_output = out_jags_mod)
}
