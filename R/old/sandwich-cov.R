#' @title Covariance Matrix of Sanwich Estimator
#' @description 
#' For incorrect likelihood function, its variance is canot be estimated using Hessian matrix,
#' thus we should use sanwich estimator to approximate the error bar.
#' @param loglik.func list of component of log-likelihood function
#' @param theta estimation of parameter
#' @export
sandwich.cov <- function(loglik.func, theta, ...) {
  p <- length(loglik.func)
  
  # define a function as the sum of loglik components
  loglik.sum <- function(theta) {
    ans <- 0
    for(ii in 1:p) {
      ans <- ans + loglik.func[[ii]](theta)
    }
    ans
  }
  
  # obtain meat part
  gradVec <- rep(NA, p)
  for(ii in 1:p) {
    gradVec[ii] <- numDeriv::grad(func = loglik.func[[ii]], x = theta)
  }
  meat <- tcrossprod(gradVec)
  
  # obtain bread part
  hess <- numDeriv::hessian(func = loglik.sum, x = theta, ...)
  bread <- solve(hess)
  
  # obtain sandwich
  bread %*% meat %*% bread
}
