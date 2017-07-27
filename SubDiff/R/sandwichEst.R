#' @title Sanwich Estimator
#' @description 
#' For incorrect likelihood function, its variance is canot be estimated using Hessian matrix, thus we 
#' use sanwich estimator to approximate the error bar.
#' @param func log-likelihood function
#' @param x estimation of parameter
#' @return vector of length of x
#' @note require package \code{numDeriv}
#' @export
sand.estimator <- function(func, x, ...) {
  p <- length(x)
  grad1 <- grad(func, x, ...)
  meat <- tcrossprod(grad1)
  hess1 <- hessian(func, x, ...)
  bread <- -1 * solve(hess1)
  var <- bread %*% meat %*% bread
  if(p == 1){
    sqrt(var[1, 1])
  } else {
    sqrt(diag(var))
  }
}