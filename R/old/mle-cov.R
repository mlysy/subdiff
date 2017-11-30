#' @title Covariance Matrix for MLE estimator
#' @description compute the covariance matrix of MLE estimator
#' @param loglik.func likelihood function
#' @param theta MLE estimation
#' @export
mle.cov <- function(loglik.func, theta, ...) {
  hess <- numDeriv::hessian(func = loglik.func, x = theta, ...)
  -1 * solve(hess)
}