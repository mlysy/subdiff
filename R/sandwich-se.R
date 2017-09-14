#' @title Sanwich Estimator
#' @description 
#' For incorrect likelihood function, its variance is canot be estimated using Hessian matrix, 
#' thus we use sanwich estimator to approximate the error bar.
#' @param func log-likelihood function
#' @param x estimation of parameter
#' @return vector of length of x
#' @note require package \code{numDeriv}
#' @export
sandwich.se <- function(func, theta, ...) {
  p <- length(theta)
  grad1 <- grad(func, theta, ...)
  meat <- tcrossprod(grad1)
  hess1 <- hessian(func, theta, ...)
  bread <- -1 * solve(hess1)
  var <- bread %*% meat %*% bread
  if(p == 1){
    ans <- sqrt(var[1, 1])
  } else {
    ans <- sqrt(diag(var))
  }
  ans
}


#' rewrite the sandwich estimator
#' its less likely to reconstructure this function as it requires uncertain number of function as 
#' input
#' Proper way of using sandwich would be write them specifically and individually.