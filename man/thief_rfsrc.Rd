% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thief_rfsrc.R
\name{thief_rfsrc}
\alias{thief_rfsrc}
\title{Temporal hierarchy reconciliation of multivariate Breiman random forests}
\usage{
thief_rfsrc(
  y,
  k = 1,
  lambda = NULL,
  frequency = 52,
  horizon = NULL,
  predict_quantiles = TRUE,
  tune_nodesize = FALSE,
  cores = parallel::detectCores() - 1,
  max_agg = NULL,
  discrete = F
)
}
\arguments{
\item{y}{\code{xts matrix}. The outcome series to be modelled. \code{NAs} are currently not supported}

\item{k}{\code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}. Default
is \code{1}, meaning that a final forecast of \code{frequency} horizons will be returned}

\item{lambda}{\code{numeric}. The Box Cox power transformation parameter for aggregate series. Must be
between \code{-1} and \code{2} inclusive. If \code{y_series} contains zeros, \code{lambda} will be set to
\code{max(c(0.7, lambda))} to ensure stability of forecasts}

\item{frequency}{\code{integer}. The seasonal frequency in \code{y}}

\item{horizon}{\code{integer}. The horizon to forecast. Defaults to \code{frequency}}

\item{predict_quantiles}{\code{logical}. If \code{TRUE}, a \code{\link[randomForestSRC]{quantreg.rfsrc}} is used
to train a second multivariate random forest using quantile loss to predict uncertainty intervals. If \code{FALSE},
the distribution of estimates from the \code{2000} original random forests is used to estimate uncertainties}

\item{tune_nodesize}{\code{logical}. If \code{TRUE}, \code{\link[randomForestSRC]{tune.nodesize}} is used to try and find
the optimal nodesize tuning parameter for the multivariate random forest. This can be slow, so the default is \code{FALSE}
and a nodesize of \code{10} is used}

\item{cores}{\code{integer}. The number of cores to use}

\item{max_agg}{(optional) \code{integer} specifying the maximum number of temporal aggregation levels
to use when reconciling, via the structural scaling method. Useful if higher levels of aggregation
are unlikely to have 'seen' recent changes in series dynamics and will likely then result in poor
forecasts as a result. Default is \code{NULL}, meaning that all levels of aggregation are used}

\item{discrete}{\code{logical} Is the series in \code{y} discrete? If \code{TRUE}, use a copula-based method
relying on the Probability Integral Transform to map the series to an approximate Gaussian distribution prior to modelling.
Forecasts are then back-transformed to the estimated discrete distribution that best fits \code{y}. Default is \code{FALSE}}
}
\value{
A \code{list} containing the reconciled forecast distributions for each series in \code{y}. Each element in
the \code{list} is a \code{horizon x 1000 matrix} of forecast predictions
}
\description{
This function fits multivariate Breiman random forests on a multivariate xts timeseries object and then
uses ensemble univariate forecast models on all higher levels of temporal aggregation to reconcile the forecasts
}
\details{
Series in \code{y} are aggregated at all possible levels up to annual using \code{\link[thief]{tsaggregates}}.
\code{\link[randomForestSRC]{rfsrc}} is used on the unaggregated series, with regressors including exponentially
weighted moving averages of each series as well as time features that can be extracted from the series. Run conditions are:
\code{ntree = 2000, nsplit = NULL, nodesize = 12}.
\code{\link{ensemble_base}} is used on all aggretated levels to find a weighted ensemble of eight
univariate forecast models that minimises mean absolute scaled error. Forecasts are then reconciled
using \code{\link[thief]{reconcilethief}} and are optionally constrained using non-negative optimisation if there
are no negative values in \code{y}. Adjustments to the original unaggregated forecast are incorporated and a
distribution of \code{1000} sample paths for each series' forecast are returned
}
\examples{
\donttest{
library(mvforecast)
data("ixodes_vets_dat")

#Fit a thief_rfsrc model
mod1 <- thief_rfsrc(y = ixodes_vets_dat$y_train,
frequency = 52, k = 1,
cores = parallel::detectCores() - 1)

#Calculate the out-of-sample CRPS
calc_crps(mod1, y_test = ixodes_vets_dat$y_test)

Plot simulation results for one of the plots in the NEON dataset
plot_mvforecast(simulation = mod1[[4]])
points(as.vector(ixodes_vets_dat$y_test[,4]))}

}
\references{
Athanasopoulos, G., Hyndman, R.,  Kourentzes, N.,  and Petropoulos, F. Forecasting with temporal hierarchies.
(2017) European Journal of Operational Research 262(1) 60–74
}
\seealso{
\code{\link{ensemble_base}}, \code{\link[randomForestSRC]{rfsrc}},
\code{\link[thief]{reconcilethief}}
}
