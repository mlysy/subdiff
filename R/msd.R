#' Sample estimate of mean squared displacement.
#'
#' @param Xt Vector or matrix of trajectories (each one is a column).
#' @param nlag Number of msd lags to calculate.
#' @param dT Optional interobservation time.  If supplied, estimates subdiffusion parameters \code{alpha, D} by regressing MSD on log-log scale.
#' @return Sample MSD vector of length \code{nlag}, or list with elements \code{msd}, \code{alpha}, \code{D}.
#' @export
msd_fit <- function(Xt, nlag, dT = NULL) {
  Xt <- as.matrix(Xt)
  N <- nrow(Xt)
  if(missing(nlag)) nlag <- floor(10 * (log10(N) - log10(ncol(Xt))))
  nlag <- min(nlag, N-1)
  # output
  ans <- rep(NA, nlag)
  for(ii in 1:nlag) {
    ans[ii] <- mean(apply(Xt, 2, diff, lag = ii)^2)
  }
  if(!is.null(dT)) {
    # log-log regression estimate
    yy <-  log(ans)
    xx <- log(1:nlag * dT)
    ybar <- mean(yy)
    xbar <- mean(xx)
    yy <- yy - ybar
    xx <- xx - xbar
    alpha <- mean(yy * xx) / mean(xx^2)
    D <- exp(ybar - alpha * xbar)
    ans <- list(msd = ans, alpha = alpha, D = D)
  }
  ans
}
