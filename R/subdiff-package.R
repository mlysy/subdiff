#' Inference for Subdiffusive Particle Tracking.
#'
#' @importFrom Rcpp evalCpp
#' @import SuperGauss
#' @import LMN
#' @importFrom numDeriv hessian
#' @importFrom stats cov optim optimize pbeta pchisq pnorm sd shapiro.test var
#' @useDynLib subdiff, .registration = TRUE
"_PACKAGE"
