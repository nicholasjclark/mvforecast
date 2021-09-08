# Estimate parameters of a Negative Binomial distribution
#'
#'This function is a wrapper for \code{\link[MASS]{fitdistr}} to generate estimates of Negative Binomial parameters
#'
#'@param x \code{vector} of values representing the discrete count vector
#'@return A \code{list} of estimates for the mean and size parameters from a \code{\link[MASS]{fitdistr}} object
#'
#'@export
#'
nb_params = function(x){
  MASS::fitdistr(x, densfun = "negative binomial")
}

# Estimate parameters of a Poisson distribution
#'
#'This function is a wrapper for \code{\link[MASS]{fitdistr}} to generate estimates of Poisson parameters
#'
#'@param x \code{vector} of values representing the discrete count vector
#'@return A \code{numeric} estimate for the mean parameter from a \code{\link[MASS]{fitdistr}} object
#'
#'@export
#'
poiss_params = function(x){
  MASS::fitdistr(x, densfun = "poisson")
}

# Transform a count vector to the Uniform distribution using the Inverse Probability Distribution Function
#'
#'For discrete series, Probability Integral Transform sampling maps observations onto a Uniform distribution that is then mapped
#'onto a Gaussian distribution for modelling. This is essentially a copula method where we assume the 'marginal'
#'is a discrete distribution. This function ranks \code{log(x + 0.01)} transformed values of \code{x} to generate
#'the Uniform(0,1) marginal and then uses \code{\link[stats]{qnorm}} to map these onto an approximately Gaussian
#'distribution
#'@param x \code{vector} of values representing the discrete count vector
#'@return A \code{vector} of approximately Gaussian PIT-transformed values
#'
#'@export
#'
paranorm = function(x){
  ranks <- data.table::frank(log2(x + 0.01), ties.method = 'max')
  pit_vals <- stats::qnorm(ranks / (length(x) + 1), sd = 2)
  return(pit_vals)
}

# Estimate copula parameters for a discrete series and return the approximately Gaussian series
#'
#'For discrete series, Probability Integral Transform sampling maps observations onto a Uniform distribution that is then mapped
#'onto a Gaussian distribution for modelling. This is essentially a copula method where we assume the 'marginal'
#'is a discrete distribution. This function determines whether a Poisson or Negative Binomial is more appropriate for the series
#'and then generates the approximately Gaussian PIT-transformed values of the series.
#'@param y \code{ts} object containing the discrete time series. \code{NA}s are allowed and will be interpolated using a combination of
#'\code{\link[zoo]{rollmean}} and \code{\link[forecast]{na.interp}}
#'@param non_neg \code{logical} indicating whether the series is restricted to be non-negative. Default is \code{TRUE}
#'@param censor \code{numeric} value ranging \code{0 - 1} indicating the upper quantile of values to truncate high outliers to
#'prior to estimating discrete distribution parameters. Useful when large outliers an lead to inflated estimates of
#'the distribution mean. Default is \code{0.99}
#'@param k \code{integer} indicating the  width of the rolling window to use for smoothly interpolating missing values. See more in
#'\code{\link[zoo]{rollmean}}
#'@return A \code{list} containing the original series with missing values interpolated, the approximately Gaussian transformed
#'series and the estimated discrete distribution parameter(s)
#'
#'@export
#'
copula_params = function(y, non_neg = TRUE, censor = 1, k = 1){
  # Remove effect of large outliers, which can lead to inflated estimates of
  # the distribution mean
  y[y > quantile(y, censor, na.rm = T)] <- quantile(y, censor, na.rm = T)

  # Use smooth interpolation of any NAs, rather than solely using forecast::na.interp (which fills NAs
  # based on STL decompositions and could be wildly inaccurate in some cases)
  y_discrete <- forecast::na.interp(zoo::rollmean(zoo::na.approx(y, na.rm = F), k = k, fill = NA))
  y_discrete <- round(as.vector(y_discrete), 0)

  if(non_neg){
    y_discrete[y_discrete < 0] <- 0
  }

  # Calculate raw discrete distribution parameters
  suppressWarnings(params <- try(nb_params(as.vector(y_discrete)), silent = TRUE))

  if(inherits(params, 'try-error')){
    suppressWarnings(params <- poiss_params(as.vector(y_discrete)))
  }

  # Get random draws from estimated distribution
  if(length(params$estimate) == 2){
    dist_mappings <- stats::qnbinom(p = seq(0, 0.99999999, length.out = 10000),
                                    size = params$estimate[1],
                                    mu = params$estimate[2])
  } else {
    dist_mappings <- stats::qpois(p = seq(0, 0.99999999, length.out = 10000),
                                  lambda = params$estimate)
  }

  # Transform y to approximately Gaussian using the PIT transformation
  para_dist <- paranorm(c(as.vector(y_discrete), dist_mappings))
  y_trans <- ts(para_dist[1:length(y_discrete)],
                start = start(y),
                frequency = frequency(y))

  return(list(y_discrete = y_discrete,
              y_trans = y_trans,
              params = params$estimate))
}

# Back-transform an approximately Gaussian series to the original discrete scale
#'
#'For discrete series, Probability Integral Transform sampling maps observations onto a Uniform distribution that is then mapped
#'onto a Gaussian distribution for modelling. This is essentially a copula method where we assume the 'marginal'
#'is a discrete distribution. This function takes an approximately Gaussian series, usually a prediction from a model,
#'and back-transforms it to the discrete distribution that was originally estimated.
#'@param x \code{vector} containing the approximately Gaussian values to be back-transformed
#'@param params the original estimated discrete distribution parameter(s). If the original was Negative Binomially distributed,
#'\code{params} will be a \code{list} of \code{length} \code{2}.
#'@return A \code{vector} containing the estimated discrete values for \code{x}
#'
#'@export
#'
back_trans = function(x, params){

  # Map forecast onto a logistic curve to respect the multiplicative nature of
  # quantile changes
  soft_boundary = function(x){
    x * pnorm(10/abs(x), sd = 2)
  }

  probs <- pnorm(as.vector(soft_boundary(x)), sd = 2)

  # No forecasts at top quantile or it will return Inf
  if(any(probs == 1)){
    probs[probs == 1] <- 0.999999
  }

   if(length(params) == 2){
          x_trans <- stats::qnbinom(p = probs,
                                    size = params[1],
                                    mu = params[2])

  } else {
        x_trans <- stats::qpois(p = probs, lambda = params[1])
  }

  if(any(is.infinite(x_trans))){
    x_trans[is.infinite(x_trans)] <- max(x_trans, na.rm = T)
  }

  return(x_trans)

}

# Back-transform an approximately Gaussian forecast to the original discrete scale
#'
#'For discrete series, Probability Integral Transform sampling maps observations onto a Uniform distribution that is then mapped
#'onto a Gaussian distribution for modelling. This is essentially a copula method where we assume the 'marginal'
#'is a discrete distribution. This function takes an approximately Gaussian forecast object
#'and back-transforms it to the discrete distribution that was originally estimated.
#'@param forecast \code{forecast} object containing the approximately Gaussian predictions to be back-transformed
#'@param orig_y \code{vector} containing the original untransformed y values
#'@param params the original estimated discrete distribution parameter(s). If the original was Negative Binomially distributed,
#'\code{params} will be a \code{list} of \code{length} \code{2}.
#'@return A \code{forecast} containing the estimated discrete predictions
#'
#'@export
#'
transform_fc_preds = function(forecast, orig_y,
                              params){
  new_fc <- forecast

  all_trans <- back_trans(c(new_fc$mean,
                            new_fc$upper[,1],
                            new_fc$upper[,2],
                            new_fc$lower[,1],
                            new_fc$lower[,2]),
                          params = params)
    indices <- c(seq(1, length(all_trans), by = length(new_fc$mean)))
    starts <- indices
    ends <- c(indices[-1]-1, length(all_trans))


    new_fc$mean <- ts(all_trans[starts[1]:ends[1]],
                      start = start(forecast$mean),
                      frequency = frequency(forecast$mean))
    new_fc$upper[,1] <- ts(all_trans[starts[2]:ends[2]],
                           start = start(forecast$mean),
                           frequency = frequency(forecast$mean))
    new_fc$upper[,2] <- ts(all_trans[starts[3]:ends[3]],
                           frequency = frequency(forecast$mean))
    new_fc$lower[,1] <- ts(all_trans[starts[4]:ends[4]],
                           start = start(forecast$mean),
                           frequency = frequency(forecast$mean))
    new_fc$lower[,2] <- ts(all_trans[starts[5]:ends[5]],
                           start = start(forecast$mean),
                           frequency = frequency(forecast$mean))
    new_fc$x <- ts(orig_y, start = start(forecast$x), frequency = frequency(forecast$x))

 return(new_fc)
}

# Calculate the highest posterior density interval
#'
#'This function uses estimated densities to calculate HPD intevals. Code originally supplied by Martyn Plummer
#'
#'@param x \code{vector} of values representing the distribution to be summarised
#'@param coverage \code{numeric} value specifying the width of the HPD interval. Default is 0.95
#'@return A \code{list} containing the reconciled forecast distributions for each series in \code{y}. Each element in
#'the \code{vector} with three values: lower estimate, median estimate and upper estimate of the HPD interval
#'
#'@export
#'
hpd <- function(x, coverage = 0.95)
{
  x <- as.matrix(x)
  out <- matrix(NA, nrow = ncol(x), ncol = 3)
  rownames(out) <- dimnames(x)[[2]]
  colnames(out) <- c("mode", "lower", "upper")

  f <- function(p) {
    if (p == density.range[2]) {
      set.coverage <- 0
    }
    else {
      p.upper <- min(y.density$y[y.density$y > p])
      p.lower <- max(y.density$y[y.density$y <= p])
      cov.upper <- sum(y.counts[y.density$y >= p.upper]) / sum(y.counts)
      cov.lower <- sum(y.counts[y.density$y >= p.lower]) / sum(y.counts)
      c <- (p.upper - p)/(p.upper - p.lower)
      set.coverage <- c * cov.upper + (1 - c) * cov.lower
    }
    return(set.coverage - coverage)
  }

  for (i in 1:ncol(x)) {
    y <- unclass(x[,i])
    y.density <- stats::density(y, n=1024)
    m <- length(y.density$x)

    ## Find the midpoint
    out[i,1] <- y.density$x[which.max(y.density$y)]
    dx <- diff(range(y.density$x)) / m
    breaks <- c(y.density$x[1] - dx / 2, y.density$x + dx /2)
    y.counts <- hist(y, breaks = breaks, plot = FALSE)$counts
    density.range <- range(y.density$y)
    uniroot.out <- stats::uniroot(f, density.range)

    ## Assuming that we have a single interval, find the limits
    out[i,2:3] <- range(y.density$x[y.density$y > uniroot.out$root])

    ## Check!
    if (sum(abs(diff(y.density$y > uniroot.out$root))) != 2) {
      warning("HPD set is not a closed interval")
    }
  }
  out <- c(out[2], out[1], out[3])
  names(out) <- c('lower', 'median', 'upper')
  return(out)
}

#' Reconcile temporal hierarchical forecasts
#'
#' A wrapper for \code{\link[thief]{reconcilethief}} that takes forecasts of time series at all levels of temporal aggregation
#' and combines them using the temporal hierarchical approach of Athanasopoulos et al (2016). Non-negative
#' optimisation constraints are used when forecasts must not be negative
#'
#' @param forecasts List of forecasts. Each element must be a time series of forecasts,
#' or a forecast object.
#' The number of forecasts should be equal to k times the seasonal period for each series,
#' where k is the same across all series.
#' @param comb Combination method of temporal hierarchies, taking one of the following values:
#' \describe{
#'   \item{"struc"}{Structural scaling - weights from temporal hierarchy}
#'   \item{"mse"}{Variance scaling - weights from in-sample MSE}
#'   \item{"ols"}{Unscaled OLS combination weights}
#'   \item{"bu"}{Bottom-up combination -- i.e., all aggregate forecasts are ignored.}
#'   \item{"shr"}{GLS using a shrinkage (to block diagonal) estimate of residuals}
#'   \item{"sam"}{GLS using sample covariance matrix of residuals}
#' }
#' @param mse A vector of one-step MSE values corresponding to each of the forecast series.
#' @param residuals List of residuals corresponding to each of the forecast models.
#' Each element must be a time series of residuals. If \code{forecast} contains a list of
#' forecast objects, then the residuals will be extracted automatically and this argument
#' is not needed. However, it will be used if not \code{NULL}.
#' @param returnall If \code{TRUE}, a list of time series corresponding to the first argument
#' is returned, but now reconciled. Otherwise, only the most disaggregated series is returned.
#' @param aggregatelist (optional) User-selected list of forecast aggregates to consider
#' @param max_agg (optional) \code{integer} specifying the maximum number of temporal aggregation levels
#' to use when reconciling, via the structural scaling method. Useful if higher levels of aggregation
#' are unlikely to have 'seen' recent changes in series dynamics and will likely then result in poor
#' forecasts as a result. Default is \code{NULL}, meaning that all levels of aggregation are used
#' @param nonnegative \code{logical} If \code{TRUE}, forecaststs are constrained using non-negative
#' optimisation to ensure negative forecasts are eliminated. Default is \code{FALSE}
#' @return
#'   List of reconciled forecasts in the same format as \code{forecast}.
#' If \code{returnall==FALSE}, only the most disaggregated series is returned.
#' @seealso \code{\link{thief}}, \code{\link{tsaggregates}}
#'
#'@export
reconcilethief_restrict <- function(forecasts,
                                  comb = c("struc","mse","ols","bu","shr","sam"),
                                  mse = NULL, residuals = NULL, returnall = TRUE,
                                  aggregatelist = NULL, max_agg = NULL,
                                  nonnegative = FALSE){

  comb <- match.arg(comb)
  if(!is.null(max_agg)){
    comb <- 'struc'
  }

  # If forecasts is a list of forecast objects, then
  # extract list of forecast time series and list of residual time series
  if(is.element("forecast",class(forecasts[[1]])))
  {
    returnclass <- "forecast"
    origf <- forecasts

    # Grab residuals
    if(is.null(residuals))
    {
      residuals <- list()
      for(i in seq_along(forecasts))
        residuals[[i]] <- residuals(forecasts[[i]])
      # Discard partial years at start of residual series
      for(i in seq_along(residuals))
      {
        tspy <- tsp(residuals[[i]])
        m <- tspy[3]
        fullyears <- trunc(length(residuals[[i]])/m)
        residuals[[i]] <- ts(utils::tail(residuals[[i]], fullyears*m),
                             frequency=m, end=tspy[2])
      }
    }

    # Grab forecasts
    for(i in seq_along(forecasts))
      forecasts[[i]] <- forecasts[[i]]$mean
  }
  else
    returnclass <- "ts"

  # Find seasonal periods
  freq <- unlist(lapply(forecasts,frequency))
  if(min(freq) != 1L)
    stop("Minimum seasonal period should be 1")
  m <- max(freq)

  # Put in order
  k <- rev(order(freq))
  forecasts <- forecasts[k]
  freq <- freq[k]

  # Check series each have same equivalent lengths.
  lengths <- unlist(lapply(forecasts, length))
  if(!is.constant(lengths/freq))
    stop("Forecast series must have the same equivalent lengths")

  # Bottom up
  if(comb == "bu")
    bts <- forecasts[[1]]

  # Some combination
  else
  {
    # Set up group matrix for hts
    # (adjusted to allow consideration of aggregatelist input)
    nsum <- rev(rep(m/freq, freq))
    unsum <- unique(nsum)
    grps <- matrix(0, nrow=length(unsum)-1, ncol=m)
    for(i in 1:(length(unsum)-1))
    {
      mi <- m/unsum[i]
      grps[i,] <- rep(1:mi, rep(unsum[i],mi))
    }

    # Set up matrix of forecasts in right structure
    nc <- length(forecasts[[1]])/m
    fmat <- matrix(0, nrow=0, ncol=nc)
    for(i in rev(seq_along(forecasts)))
      fmat <- rbind(fmat, matrix(forecasts[[i]], ncol=nc))

    # OLS or WLS reconciliation
    if(is.element(comb, c("struc","ols","mse")))
    {
      if(comb=="struc")
        weights <- 1/nsum
      else if(comb=="ols")
        weights <- NULL
      else if(comb=="mse")
        weights <- 1/rep(rev(mse), rev(unsum))

      # Restrict levels of aggregation to use based on max_agg by setting levels above
      # max_agg as zero weights in the structural method
      if(!is.null(max_agg)){
        tot_levels_agg <- length(unique(nsum))
        if(max_agg > tot_levels_agg){
          warning('Max agg greater than maximum level of aggregation. Ignoring')
          max_agg <- max(nsum)
        } else {
          max_agg <- rev(unique(nsum))[max_agg]
        }
        weights[nsum > max_agg] <- 0
      }

      bts <- hts::combinef(t(fmat), groups=grps, weights=weights, keep='bottom',
                           nonnegative = nonnegative)
    }
    else
      # GLS reconciliation
    {
      # Set up matrix of residuals in right structure
      if(is.null(residuals))
        stop("GLS needs residuals")
      nc <- length(residuals[[1]])/m
      rmat <- matrix(0, nrow=0, ncol=nc)
      for(i in rev(seq_along(forecasts)))
        rmat <- rbind(rmat, matrix(residuals[[i]], ncol=nc))
      bts <- hts::MinT(t(fmat), groups=grps, residual=t(rmat),
                       covariance=comb, keep='bottom', nonnegative = nonnegative)
    }
    # Turn resulting reconciled forecasts back into a ts object
    bts <- ts(c(t(bts)))
    tsp(bts) <- tsp(forecasts[[1]])
  }

  # Now figure out what to return
  if(returnclass == "ts") # Just return time series
  {
    if(!returnall)
      return(bts)
    else
      return(thief::tsaggregates(bts, aggregatelist = aggregatelist))
  }
  else
    # Return forecast objects
  {
    if(!returnall)
    {
      adj <- bts - origf[[1]]$mean
      origf[[1]]$mean <- bts
      if(!is.null(origf[[1]]$lower))
      {
        origf[[1]]$lower <- sweep(origf[[1]]$lower,1,adj,"+")
        origf[[1]]$upper <- sweep(origf[[1]]$upper,1,adj,"+")
      }
      return(origf[[1]])
    }
    else
    {
      allts <- thief::tsaggregates(bts, aggregatelist=aggregatelist)
      for(i in seq_along(origf))
      {
        adj <- allts[[i]] - origf[[i]]$mean
        origf[[i]]$mean <- allts[[i]]
        if(!is.null(origf[[i]]$lower))
        {
          origf[[i]]$lower <- sweep(origf[[i]]$lower,1,adj,"+")
          origf[[i]]$upper <- sweep(origf[[i]]$upper,1,adj,"+")
        }
      }
      return(origf)
    }
  }
}
