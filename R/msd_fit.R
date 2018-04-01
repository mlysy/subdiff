#' Sample estimate of mean squared displacement.
#'
#' @param Xt Vector or matrix of trajectories (each one is a column).
#' @param nlag Number of MSD lags to calculate.
#' @param demean Logical; whether or not to remove the mean of \code{Xt} as estimated by linear drift.
#' @param dT Optional interobservation time.  If supplied, estimates subdiffusion parameters \code{(alpha, D)} by regressing MSD over time on log-log scale.
#' @param msd Optionally pre-computed sample MSD, in which case it is used for regression estimate of \code{(alpha, D)}.
#' @return Sample MSD vector of length \code{nlag}, or list with elements \code{msd}, \code{alpha}, \code{D}.
#' @export
msd_fit <- function(Xt, nlag, demean = TRUE, dT = NULL, msd) {
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
    msd <- rep(NA, nlag)
    for(ii in 1:nlag) {
      msd[ii] <- mean(apply(Xt, 2, diff, lag = ii)^2)
    }
  } else {
    nlag <- length(msd)
  }
  if(!is.null(dT)) {
    # log-log regression estimate
    yy <-  log(msd)
    xx <- log(1:nlag * dT)
    ybar <- mean(yy)
    xbar <- mean(xx)
    yy <- yy - ybar
    xx <- xx - xbar
    alpha <- mean(yy * xx) / mean(xx^2)
    D <- exp(ybar - alpha * xbar)
    msd <- list(msd = msd, alpha = alpha, D = D)
  }
  msd
}
