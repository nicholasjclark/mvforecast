#'Fit a Bayesian forecast model to a multivariate set of discrete time series
#'
#'This function takes fitted values from a log-transformed forecast models and uses those as
#'predictors in a Negative Binomial regression to model the discrete outcome time series. Multivariate
#'dependencies are incorporated via a covariance matrix
#'
#'
#'@param y \code{xts matrix}. The discrete outcome series to be modelled. \code{NAs} are allowed
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
#'@details The discrete series are first log-transformed and interpolated to fit univariate forecast models.
#'Fitted values from the forecast models are then used as predictors in a Negative Binomial Bayesian regression
#'to model the raw discrete series, while inter-series dependencies are captured using a multivariate
#'covariance matrix
#'
#'@export
#'
fc_jags_mv = function(y,
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
  # Variable checks
  if(use_nb){
    use_nb = 1
  } else {
     use_nb = 0
  }

  n <- NCOL(y)
  if (n == 1){
    stop("\ncannot specify a multivariate model with only one series")
  }

  ynames <- colnames(y)
  if(is.null(ynames)) {
    colnames(y) <- paste0("Series", 1:n)
    ynames <- colnames(y)
  }


  if(!tolower(model) %in% c('tbats', 'arima', 'stlf')){
    stop("model must be either 'TBATS', 'arima' or 'stlf'")
  }

  # Set forecast horizon if missing
  if(missing(horizon)){
    horizon <- frequency
  }

  # Generate univariate forecast model fitted values
  fc_fits <- do.call(cbind, lapply(seq_len(n), function(series){
    if(tolower(model) == 'tbats'){
      fc_fit <- fit_tbats_discreet(y[,series],
                                   frequency = frequency,
                                   horizon = horizon)
    }

    if(tolower(model) == 'arima'){
      fc_fit <- fit_arima_discreet(y[,series],
                                   frequency = frequency,
                                   horizon = horizon)
    }

    if(tolower(model) == 'stlf'){
      fc_fit <- fit_stlf_discreet(y[,series],
                                  frequency = frequency,
                                  horizon = horizon)
    }
    fc_fit
  }))


y_orig <- y
y <- rbind(as.matrix(y), matrix(NA, ncol = NCOL(y), nrow = horizon))

model_file <- " model {

## Mean expectations
 for (i in 1:n) {
  for (s in 1:n_series){
   mu[i, s] <- alpha[s] + fc_fits[i, s];
  }
 }

for (i in 1:n) {
 Z[i, 1:n_series] ~ dmnorm(mu[i, 1:n_series], tau[ , ])
}

for (s in 1:n_series){
 alpha_raw[s] ~ dnorm(0, 1)
 alpha[s] <- mu_alpha + (sigma_alpha * alpha_raw[s])
}
# alpha should be negative as we added 1/6 prior to the univariate model
mu_alpha ~ dunif(-2, 0)
sigma_alpha ~ dgamma(0.1, 0.1)

## Negative binomial likelihood
for (i in 1:n) {
 for (s in 1:n_series){
  y[i, s] ~ dnegbin(rate[i, s], r);
  rate[i, s] <- ifelse((r / (r + exp(Z[i, s]))) < min_eps, min_eps,
                       (r / (r + exp(Z[i, s]))))
 }
}

## Set overdispersion parameter to 10000 if not using nb;
## this approximates the Poisson distribution
 r1 ~ dgamma(0.01, 0.01)
 r2 <- 10000
 r_indicator <- ifelse(use_nb > 0, 1, 0)
 r <- (r1 * r_indicator) + (r2 * (1 - r_indicator))

## Wishart is appropriate conjugate prior for dmnorm precision
# df = n_series + 1 sets uniform distribution on correlation parameters
tau[1:n_series, 1:n_series] ~ dwish(diag_mat[ , ], n_series + 1)

# Scaled residual correlation matrix (between -1 and 1; Gelman & Hill 2007)
Covar.raw[1:n_series, 1:n_series] <- inverse(tau[ , ])
 for(i in 1:n_series){
  for(j in 1:n_series){
   cor[i, j] <- Covar.raw[i, j] / sqrt(Covar.raw[i, i] * Covar.raw[j, j])
  }
}

 ## Posterior predictions
 for (i in 1:n) {
  for (s in 1:n_series){
   ypred[i, s] ~ dnegbin(rate[i, s], r)
  }
}

}"

jags_data <- list(
  y = y,
  n_series = NCOL(y),
  fc_fits = fc_fits,
  n = NROW(y),
  diag_mat = diag(NCOL(y)),
  min_eps = .Machine$double.eps,
  use_nb = use_nb
)

inits <- function(){
  list(alpha_raw = rep(-0.5, NCOL(y)))
}

# Initiate adaptation of the model
load.module("glm")
jags_mod <- jags.model(textConnection(model_file),
                      data = jags_data,
                      inits = inits,
                      n.chains = n.chains,
                      n.adapt = n.adapt)


# Update the model for the burnin period
update(jags_mod, n.burnin)

# Sample from the posterior
out_jags_mod <- rjags::coda.samples(jags_mod,
                              variable.names = c('alpha',
                                                 'ypred',
                                                 'r',
                                                 'cor'),
                              n.iter = n.iter,
                              thin = thin)

if(auto_update){
  # Update until reasonable convergence in the form of Rhat and ESS
    update_params <- c('alpha','r')

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
                                                          'r', 'cor'),
                                       n.iter = n.iter,
                                       thin = thin)
    mod_summary <- MCMCvis::MCMCsummary(out_jags_mod, update_params)
  }
}

# Return forecast and fitted jags model
ends <- seq(0, dim(MCMCvis::MCMCchains(out_jags_mod, 'ypred'))[2],
                    length.out = NCOL(y) + 1)
starts <- ends + 1
starts <- c(1, starts[-c(1, (NCOL(y)+1))])
ends <- ends[-1]

fc_list <- lapply(seq_len(n), function(series){
  ypreds <- MCMCvis::MCMCchains(out_jags_mod, 'ypred')[,starts[series]:ends[series]]
  t(ypreds[, (NROW(y_orig)+1):NCOL(ypreds)])
})
names(fc_list) <- ynames


 list(forecast = fc_list, jags_output = out_jags_mod)
}
