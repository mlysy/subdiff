#' @title Wrapper function for Mode polynomial
#' @description 
#' MODEPOLY given a polynomial with distinct roots a_0 < ... < a_N, find its
#' local modes using the golden section method programmed in c++.
#' @param lambda
#' @param nSteps number of steps
#' @param tol tolerence 
#' @return vector of length of \code{lambda}
#' @export
modePoly <- function(lambda, nSteps, tol){
  ModePoly(lambda, nSteps, tol)
}