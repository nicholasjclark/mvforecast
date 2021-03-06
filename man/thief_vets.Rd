% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thief_vets.R
\name{thief_vets}
\alias{thief_vets}
\title{Temporal hierarchy reconciliation of vector exponential smoothing models}
\usage{
thief_vets(
  y,
  multi_freq = 12,
  k = 1,
  level = "grouped",
  slope = "none",
  damped = "none",
  seasonal = "common",
  lambda = NULL,
  dependence = "equicorrelation",
  frequency = 52,
  horizon = NULL,
  cores = parallel::detectCores() - 1,
  group = NULL,
  xreg = NULL,
  max_agg = NULL,
  xreg_include = NULL,
  newxreg = NULL,
  save_plots = FALSE,
  fig_path = "",
  discrete = FALSE
)
}
\arguments{
\item{y}{\code{xts matrix}. The outcome series to be modelled. \code{NAs} are currently not supported}

\item{multi_freq}{\code{integer}. Minimum frequency in the temporal hierarchy that will be modelled using
a vector exponential smoothing model. Aggregates with frequencies below this threshold will be modelled with
univariate models using the \code{\link[forecast]{forecast}} function}

\item{k}{\code{integer} specifying the length of the forecast horizon in multiples of \code{frequency}. Default
is \code{1}, meaning that a final forecast of \code{frequency} horizons will be returned}

\item{level}{\code{character}. The dynamics for the level component in the \code{tsvets} model. Options are:
c("constant", "diagonal", "common", "full", "grouped")}

\item{slope}{\code{character}. The dynamics for the slope component in the \code{tsvets} model. Options are:
c("none", "constant", "common", "diagonal", "full", "grouped")}

\item{damped}{\code{character}. The dynamics for the dampening component in the \code{tsvets} model. Options are:
c("none", "common", "diagonal", "full", "grouped")}

\item{seasonal}{\code{character}. The dynamics for the seasonal component in the \code{tsvets} model. Options are:
c("none", "common", "diagonal", "full", "grouped")}

\item{lambda}{\code{numeric proportional}. The multivariate Box Cox power transformation parameter for all series.
Must be between \code{0} and \code{1.5} inclusive}

\item{dependence}{\code{character}. The multivariate error dependence structure to impose. Options are:
c("diagonal", "full", "equicorrelation", "shrinkage")}

\item{frequency}{\code{integer}. The seasonal frequency in \code{y}}

\item{horizon}{\code{integer}. The horizon to forecast. Defaults to \code{frequency}}

\item{cores}{\code{integer}. The number of cores to use. This is used to initialize the states of each series
using \code{\link[tsets]{ets_modelspec}}}

\item{group}{Optional \code{vector} of indices denoting which group the series belongs to (when using the grouped dynamics).
Defaults to \code{NULL}}

\item{xreg}{Optional \code{xts} matrix of external regressors. Defaults to \code{NULL}}

\item{max_agg}{(optional) \code{integer} specifying the maximum number of temporal aggregation levels
to use when reconciling, via the structural scaling method. Useful if higher levels of aggregation
are unlikely to have 'seen' recent changes in series dynamics and will likely then result in poor
forecasts as a result. Default is \code{NULL}, meaning that all levels of aggregation are used}

\item{xreg_include}{optional \code{matrix} of dimension \code{ncol(y)} by \code{ncol(xreg)} populated
with either 0, 1 or 2+ (0 = no beta, 1 = individual beta and 2 = grouped beta). It is also
possible to have group wise pooling. For instance 2 variables sharing one pooled estimates,
and 3 other variables sharing another grouped estimate would have values of (2,2,3,3,3).
The index for group wise pooling starts at 2 and should be incremented for each new group added. Defaults to \code{NULL}}

\item{newxreg}{Optional \code{xts} matrix of future values for external regressors to be used in forecasting.
Defaults to \code{NULL}}

\item{save_plots}{Logical. Plots of fitted and residual values for the unaggregated series will be saved to
\code{fig_path} if \code{TRUE}}

\item{fig_path}{\code{character}. Optional filepath where fitted and residual plots will be saved. Defaults to the
current working directory}

\item{discrete}{\code{logical} Is the series in \code{y} discrete? If \code{TRUE}, use a copula-based method
relying on the Probability Integral Transform to map the series to an approximate Gaussian distribution prior to modelling.
Forecasts are then back-transformed to the estimated discrete distribution that best fits \code{y}. Default is \code{FALSE}}
}
\value{
A \code{list} containing the reconciled forecast distributions for each series in \code{y}. Each element in
the \code{list} is a \code{horizon x 1000 matrix} of forecast predictions
}
\description{
This function fits vector exponential smoothing on the base level and optionally on
higher levels of temporal aggregation for a multivariate xts timeseries object
}
\details{
Series in \code{y} are aggregated at all possible levels up to annual using \code{\link[thief]{tsaggregates}}.
\code{\link[tsvets]{estimate.tsvets.spec}} is used on the unaggregated (original series), and optionally on higher levels
of aggregation down to a frequency of \code{multi_freq}. At frequencies below \code{multi_freq}, the best-fitting
univariate model is chosen using automatic ensembles in \code{\link{ensemble_base}}. Forecasts are reconciled
using \code{\link[thief]{reconcilethief}} and are optionally constrained using non-negative optimisation if \code{lambda}
is provided. Adjustments to the original unaggregated forecast are incorporated and a distribution of \code{1000} sample
paths for each series' forecast are returned
}
\examples{
\donttest{
library(mvforecast)
data("ixodes_vets_dat")

# View the returned data
head(ixodes_vets_dat$y_train)
head(ixodes_vets_dat$xreg_train)
ixodes_vets_dat$xreg_include
head(ixodes_vets_dat$future_xreg)

# Fit a vets model with no regressors and common seasonality with the tsvets package
mod1 <- tsvets:::simulate.tsvets.estimate(tsvets:::estimate.tsvets.spec(tsvets:::vets_modelspec(ixodes_vets_dat$y_train,
level = "grouped",
slope = "none",
damped = "none",
seasonal = "common",
lambda = 1,
dependence = "equicorrelation",
frequency = 52,
cores = parallel::detectCores() - 1,
group = ixodes_vets_dat$groups),
solver = "solnp",
control = list(trace = 0)),
nsim = 1000,
h = ixodes_vets_dat$h)

# Calculate a summary of the summed out-of-sample CRPS
calc_crps(simulation = mod1, y_test = ixodes_vets_dat$y_test)

# Explore whether reconciliation of temporal hierarchies improves predictions
mod2 <- thief_vets(y = ixodes_vets_dat$y_train,
multi_freq = 12,
level = "grouped",
slope = "none",
damped = "none",
seasonal = "common",
lambda = 1,
dependence = "equicorrelation",
frequency = 52,
cores = parallel::detectCores() - 1,
group = ixodes_vets_dat$groups,
save_plots = FALSE)

calc_crps(simulation = mod2, y_test = ixodes_vets_dat$y_test)

# Plot one of the forecasts against the true values in the test set
plot_mvforecast(simulation = mod2[[4]])
points(as.vector(ixodes_vets_dat$y_test[,4]))}

}
\references{
Athanasopoulos, G and de Silva, A. (2012),Multivariate Exponential Smoothing for Forecasting Tourist Arrivals,
Journal of Travel Research 51(5) 640–652.\cr\cr
de Silva, A., R. Hyndman, and R. D. Snyder. (2010).The Vector Innovations Structural Time Series Framework: A Simple Approach
to Multivariate Forecasting, Statistical Modelling (10) 353–74.\cr\cr
Athanasopoulos, G., Hyndman, R.,  Kourentzes, N.,  and Petropoulos, F. Forecasting with temporal hierarchies.
(2017) European Journal of Operational Research 262(1) 60–74
}
\seealso{
\code{\link[tsvets]{vets_modelspec}}, \code{\link[tsets]{ets_modelspec}}, \code{\link{ensemble_base}},
\code{\link[thief]{reconcilethief}}
}
