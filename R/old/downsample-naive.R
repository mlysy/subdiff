#' @title fBM inference using composite downsampling
#' @description We assume that time series dY follows MN(Beta * dT, V_incr.fBM, Sigma),
#' thus downsampled dY_ds follows MN(Beta * ds * dT, V_ds.incr.fBM, Sigma).
#' @param Xt Time series for inference
#' @param dT Interobservation time
#' @param Tz Toeplitz-object of size floor(nrow(Xt)/ds)-1
#' @param ds Downsampling rate
#' @param pos position of downsampling, integer between 1 and ds, default to be 1.
#' @param varCal logic, return the covariance matrix of MLE estimate if true
#' @param se logic, return the standard error of MLE instead of covariance matix if true
#' @param trans logic, return the 'normalized' parameters
#' @export
downsample.naive <- function(Xt, dT, Tz, ds, pos = 1, varCal = FALSE, se = FALSE, trans = FALSE) {
  Xt <- downSample(Xt, ds, pos)
  dX <- apply(Xt, 2, diff)
  N <- nrow(dX)
  q <- ncol(dX)
  dT <- dT*ds

  if(missing(Tz)) {
    Tz <- Toeplitz(n = N)
  }else if(nrow(Tz) != N) {
    stop("Tz has incompatible dimension with Xt and ds")
  }

  alpha <- optimize(f = function(alpha) {
    Tz$setAcf(fbm_acf(alpha, dT, N))
    lmn.prof(Y = dX, X = dT, acf = Tz)
  }, interval = c(0, 2), maximum = TRUE)$maximum

  acf1 <- fbm_acf(alpha, dT, N)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  ans <- theta <- transFunc(alpha, suff$Beta.hat, suff$S/N, trans)
  if(varCal) {
    vcov <- mle.cov(loglik.func = function(theta) {
      suff <- itransFunc(theta, trans)
      Tz$setAcf(fbm_acf(suff$alpha, dT = dT, N = N))
      composite.full(Y = Xt, X = dT, Beta = suff$Beta, Sigma = suff$Sigma, acf = Tz, ds = 1)
    }, theta = theta)
    if(se) {
      vcov <- sqrt(diag(vcov))
    }
    ans <- list(mle = theta, cov = vcov)
  }
  ans
}
