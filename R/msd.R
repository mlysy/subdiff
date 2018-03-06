#' @title Mean Square Displacement
#' @param Xt trajectory
#' @param dT interobservation time, required when computing slope \code{alpha}
#' @param nlag order of MSD of interested, default 20
#' @param msd.only logical, return both computed msd and slope of msd when TRUE
#' @return vector of length \code{nlag}
#' @export
msd <- function(Xt, dT, nlag = 30, msd.only = TRUE) {
  N <- length(Xt)
  ans <- rep(NA, nlag)
  for(ii in 1:nlag) {
    ans[ii] <- mean(diff(Xt, lag = ii)^2)
  }
  if(!msd.only) {
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
