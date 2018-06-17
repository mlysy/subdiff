#' Sample estimate of mean squared displacement.
#'
#' @param Xt Vector or matrix of trajectories (each one is a column).
#' @param nlag Number of MSD lags to calculate.
#' @param demean Logical; whether or not to remove the mean of \code{Xt} as estimated by linear drift.
#' @param dT Optional interobservation time.  If supplied, estimates subdiffusion parameters \code{(alpha, D)} by regressing MSD over time on log-log scale.
#' @param logw Logical; whether or not regression should be log-weighted (see Details).
#' @param msd Optionally pre-computed sample MSD, in which case it is used for the regression estimate of \code{(alpha, D)}.
#' @return Sample MSD vector of length \code{nlag}, or list with elements \code{msd}, \code{alpha}, \code{D}.
#' @details When \code{dT} is supplied, a power-law fit to the empirical MSD is calculated by performing the linear regression
#' \preformatted{
#' log(msd) ~ log(D) + alpha * log(time).
#' }
#' From a curve fitting perspective, we seek to estimate the "best" power-law describing the MSD, i.e., if MSD(t) is the true MSD at time t, then the best power law minimizes
#' \preformatted{
#' integral_0^Inf | log MSD(t) - log(D) - alpha * log(t) |^2 w(t) dt,
#' }
#' where \code{w(t) > 0} is a weight function.  The default value of \code{logw = FALSE} uses \code{w(t) = 1}, but this doesn't give the best MSD fit on the log-log scale, as there are exponentially more points as we move right in the graph, such that the right side of the graph will dominate the fit.  Setting \code{logw = TRUE} uses uniform weighting on the log-scale, which corresponds to \code{w(t) = 1/t}.
#' @export
msd_fit <- function(Xt, nlag, demean = TRUE, dT = NULL, logw = FALSE, msd) {
  if(missing(msd)) {
    Xt <- as.matrix(Xt)
    N <- nrow(Xt)
    if(demean) {
      mu <- (Xt[N,] - Xt[1,])/(N-1)
      Xt <- Xt - ((1:N-1) %o% mu)
    }
    if(missing(nlag)) nlag <- floor(10 * (log10(N) - log10(ncol(Xt))))
    nlag <- min(nlag, N-1)
    # output
    msd <- .SampleMSD(Xt, nlag)
    ## msd <- rep(NA, nlag)
    ## for(ii in 1:nlag) {
    ##   msd[ii] <- mean(apply(Xt, 2, diff, lag = ii)^2)
    ## }
  } else {
    nlag <- length(msd)
  }
  if(!is.null(dT)) {
    ww <- if(logw) 1/(1:nlag) else rep(1, nlag)
    ww <- ww/sum(ww)
    # log-log regression estimate
    yy <-  log(msd)
    xx <- log(1:nlag * dT)
    ybar <- sum(ww * yy)
    xbar <- sum(ww * xx)
    yy <- yy - ybar
    xx <- xx - xbar
    alpha <- sum(ww * yy * xx) / sum(ww * xx^2)
    D <- exp(ybar - alpha * xbar)
    msd <- list(msd = msd, alpha = alpha, D = D)
  }
  msd
}
