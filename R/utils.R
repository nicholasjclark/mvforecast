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
#'
#' @return
#'   List of reconciled forecasts in the same format as \code{forecast}.
#' If \code{returnall==FALSE}, only the most disaggregated series is returned.
#' @seealso \code{\link{thief}}, \code{\link{tsaggregates}}
#'
#'@export
reconcilethief_nonneg <- function(forecasts,
                                  comb = c("struc","mse","ols","bu","shr","sam"),
                                  mse = NULL, residuals = NULL, returnall = TRUE,
                                  aggregatelist = NULL){

  comb <- match.arg(comb)

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
      bts <- hts::combinef(t(fmat), groups=grps, weights=weights, keep='bottom', nonnegative = T)
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
                       covariance=comb, keep='bottom', nonnegative = T)
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
