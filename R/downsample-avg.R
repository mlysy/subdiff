#' @title fBM inference using average downsampling estimator
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
downsample.avg <- function(Xt, dT, Tz, ds, varCal = FALSE, se = FALSE, trans = FALSE) {
  N <- nrow(Xt)
  q <- ncol(Xt)
  N.ds <- floor(N/ds) - 1

  if(q == 1) {
    thetaMat <- matrix(NA, ds, 3)
    vcov <- matrix(NA, 3, 3)
  } else if (q == 2) {
    thetaMat <- matrix(NA, ds, 6)
    vcov <- matrix(0, 6, 6)
  } else {
    stop("So far we can only deal with dimension < 3")
  }

  if(missing(Tz)) {
    Tz <- Toeplitz(n = N.ds)
  }else if(nrow(Tz) != N.ds) {
    stop("Tz has incompatible dimension with Xt and ds")
  }


  for(ii in 1:ds) {
    thetaMat[ii, ] <- downsample.naive(Xt, dT, Tz, ds, pos = ii, trans = trans)
  }

  ans <- theta <- apply(thetaMat, 2, mean)

  if(varCal) {
    for(ii in 1:ds) {
      Xt.ds <- downSample(Xt, ds, pos = ii)
      dX <- apply(Xt.ds, 2, diff)
      vcov <- vcov + mle.cov(loglik.func = function(theta) {
        suff <- itransFunc(theta, trans)
        Tz$setAcf(fbm_acf(suff$alpha, dT = dT, N = N.ds))
        composite.full(Y = Xt.ds, X = ds*dT, Beta = suff$Beta, Sigma = suff$Sigma, acf = Tz, ds = 1)
      }, theta = thetaMat[ii, ])
    }
    vcov <- vcov / ds^2
    if(se) {
      vcov <- sqrt(diag(vcov))
    }
    ans <- list(mle = theta, cov = vcov)
  }
  ans
}
