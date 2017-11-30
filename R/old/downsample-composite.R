#' @title fBM inference using composite downsampling
#' @description We assume that time series dY follows MN(Beta * dT, V_incr.fBM, Sigma),
#' thus downsampled dY_ds follows MN(Beta * ds * dT, V_ds.incr.fBM, Sigma).
#' @param Xt Time series for inference
#' @param dT Interobservation time
#' @param Tz Toeplitz-object of size floor(nrow(Xt)/ds)-1
#' @param ds Downsampling rate
#' @param varCal logic, return the covariance matrix of MLE estimate if true
#' @param se logic, return the standard error of MLE instead of covariance matix if true
#' @param trans logic, return the 'normalized' parameters
#' @export
downsample.composite <- function(Xt, dT, Tz, ds, varCal = FALSE, se = FALSE, trans = FALSE) {
  N <- nrow(Xt)
  q <- ncol(Xt)
  N.ds <- floor(N/ds) - 1

  if(missing(Tz)) {
    Tz <- Toeplitz(n = N.ds)
  }else if(nrow(Tz) != N.ds) {
    stop("Tz has incompatible dimension with Xt and ds")
  }

  alpha <- optimize(f = function(alpha) {
    Tz$setAcf(fbm_acf(alpha, dT*ds, floor(nrow(Xt)/ds) - 1))
    composite.prof(Y = Xt, X = dT, acf = Tz, ds = ds)
  }, interval = c(0, 2), maximum = TRUE)$maximum

  acf1 <- fbm_acf(alpha, dT*ds, N.ds)
  Tz$setAcf(acf1)
  suff <- composite.suff(Y = Xt, X = dT, acf = Tz, ds = ds)
  ans <- theta <- transFunc(alpha, suff$Betahat, suff$S, trans)
  if(varCal) {
    vcov <- mle.cov(loglik.func = function(theta) {
      suff <- itransFunc(theta, trans)
      Tz$setAcf(fbm_acf(suff$alpha, dT = dT*ds, N = floor(nrow(Xt)/ds) - 1))
      composite.full(Y = Xt, X = dT, Beta = suff$Beta, Sigma = suff$Sigma, acf = Tz, ds = ds)
    }, theta = theta)
    if(se) {
      vcov <- sqrt(diag(vcov))
    }
    ans <- list(mle = theta, cov = vcov)
  }
  ans
}

#' nested function for profile composite likelihood
# .fbm.prof.cl <- function(alpha, Xt, dT, ds, Tz) {
#   N <- floor(nrow(Xt)/ds) - 1
#   acf1 <- fbm_acf(alpha, dT * ds, N)
#   Tz$setAcf(acf1)
#   composite.prof(Y = Xt, X = dT, acf = Tz, ds = ds)
# }
#
# fbm.full.ll <- function(theta, dX, dT, ds, Tz, trans) {
#   N <- nrow(dX)
#   ans <- itransFunc(theta, trans)
#   alpha <- ans$alpha
#   Sigma <- ans$Sigma
#   Beta <- ans$Beta
#   dG <- matrix(dT, N, 1)
#   dX <- dX - dG %*% Beta
#   acf1 <- fbm_acf(alpha, dT, N)
#   Tz$setAcf(acf1)
#   ldV <- determinant(Tz)
#   ldS <- log(det(Sigma))
#   ll <- N * ldS + 2 * ldV + tr(solve(Sigma) %*% crossprod(dX, solve(Tz, dX)))
#   -ll / 2
# }
