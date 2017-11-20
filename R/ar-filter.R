#' @title AR(1) Filter
#' @description Model: dY(t) = (1-rho) dX(t) + rho dY(t-1)
#' @param dX Time series for inference
#' @param dT Interobservation time
#' @param Tz Toeplitz-object of size N
#' @param varCal logic, return the covariance matrix of MLE estimate if true
#' @param se logic, return the standard error of MLE instead of covariance matix if true
#' @param trans logic, return the 'normalized' parameters
#' @export
ar.filter <- function(dX, dT, Tz, varCal = FALSE, se = FALSE, trans = FALSE) {
  N <- nrow(dX)
  q <- ncol(dX)

  if(missing(Tz)) {
    Tz <- Toeplitz(n = N)
  }

  theta <- optim(par = c(1, 0), fn = function(theta) {
    alpha <- theta[1]
    rho <- theta[2]
    if(rho <= -1 || rho >= 1) {
      -Inf
    } else {
      dY <- ar2eps(dX, rho)
      Tz$setAcf(fbm_acf(alpha, dT, N))
      lmn.prof(Y = dY, X = dT, acf = Tz) - log(1-rho)*N*q
    }
  }, control = list(fnscale = -1))$par

  Tz$setAcf(fbm_acf(theta[1], dT, N))
  dY <- ar2eps(dX, theta[2])
  suff <- lmn.suff(Y = dY, X = dT, acf = Tz)

  # order of theta: 1 alpha, 2:3 beta, 4:6 sigma, 7 rho
  ans <- theta <- c(transFunc(theta[1], suff$Beta.hat, suff$S/N, trans), theta[2])

  if(varCal) {
    vcov <- mle.cov(loglik.func = function(theta) {
      suff <- itransFunc(theta[1:6], trans)
      Tz$setAcf(fbm_acf(suff$alpha, dT = dT, N = N))
      dY <- ar2eps(dX, theta[7])
      Yt <- apply(rbind(0, dY), 2, cumsum)
      composite.full(Y = Yt, X = dT, Beta = suff$Beta, Sigma = suff$Sigma, acf = Tz, ds = 1)
    }, theta = theta)
    if(se) {
      vcov <- sqrt(diag(vcov))
    }
    ans <- list(mle = theta, cov = vcov)
  }
  ans
}

#' nested function
#' can we export Rcpp function with roxgen2?
#' convert AR(1) process into its eps with given rho
ar2eps <- function(Yt, rho){
  N <- nrow(Yt)
  Xt <- Yt / (1-rho)
  Xt[-1, ] <- (Yt[-1, ] - rho * Yt[-N, ]) / (1-rho)
  Xt
}
