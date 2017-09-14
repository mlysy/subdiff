#' @title Standard Error for Normal Likelihood
#' @description 
#' @details 
#' @param func likelihood function
#' @param theta MLE estimation
#' @export
loglik.se <- function(func, theta, ...) {
  p <- length(theta)
  hess1 <- hessian(func = func, x = theta, ...)
  var <- -1 * solve(hess1)
  if(p == 1) {
    ans <- sqrt(var)
  } else {
    ans <- sqrt(diag(var))
  }
  ans
}