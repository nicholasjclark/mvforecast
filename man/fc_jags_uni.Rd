% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_jags_uni.R
\name{fc_jags_uni}
\alias{fc_jags_uni}
\title{Fit a Bayesian forecast model to a univariate discrete time series}
\usage{
fc_jags_uni(
  y,
  model = "arima",
  horizon,
  frequency,
  use_nb = TRUE,
  n.chains = 2,
  n.adapt = 1000,
  n.burnin = 10000,
  n.iter = 10000,
  thin = 10,
  auto_update = TRUE
)
}
\arguments{
\item{y}{\code{xts vector}. The discrete outcome series to be modelled. \code{NAs} are allowed}

\item{model}{\code{character} either 'TBATS', 'arima' or 'stlf'}

\item{horizon}{\code{integer}. The horizon to forecast. Defaults to \code{frequency}}

\item{frequency}{\code{integer}. The seasonal frequency in \code{y}}

\item{use_nb}{\code{logical} If \code{TRUR}, use a Negative Binomial likelihood with estimated
overdispersion parameter r;}

\item{n.chains}{\code{integer} the number of parallel chains for the model}

\item{n.adapt}{\code{integer} the number of iterations for adaptation. See adapt for details.
If \code{n.adapt} = 0 then no adaptation takes place}

\item{n.iter}{\code{integer} the number of iterations of the Markov chain to run}

\item{auto_update}{\code{logical} If \code{TRUE}, the model is run for up to \code{3} additional sets of
\code{n.iter} iterations, or until the lower 15th percentile of effective sample sizes reaches \code{200}}
}
\value{
A \code{list} containing the posterior forecast from the forecast model and the original fitted
\code{\link[rjags]{jags.model}}
}
\description{
This function takes fitted values from a log-transformed forecast model and uses those as
a predictor in a Negative Binomial regression to predict the discrete outcome time series
}
\details{
The discrete series is first log-transformed and interpolated to fit a forecast model. Fitted values
from the forecast model are then used as a predictor in a Negative Binomial Bayesian regression to model
the raw discrete series
}
