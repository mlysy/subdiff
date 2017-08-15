#' @title Fisher Information Estimation for MLE variance
#' @description MLE estimation asymptotically follows normal distribution with mean unbiased and 
#' variance Fisher Information, which equals the negative inverse of Hessian matrix under 
#' most circumstance.
#' @param func log-likelihood function
#' @param x estimation of parameter
#' @return vector of length of x
#' @note require package \code{numDeriv}
#' @export
hess.estimator <- function(func, x, ...) {
  p <- length(x)
  hess1 <- hessian(func, x, ...)
  var <- -1 * solve(hess1)
  if(p == 1) {
    sqrt(var[1, 1])
  } else {
    sqrt(diag(var))
  }
}