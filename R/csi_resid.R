#' Residuals of a Gaussian CSI process.
#'
#' @template args-dX
#' @template args-dT
#' @param mu Vector of linear drift in each direction.
#' @param acf Autocorrelation, a vector of length `nrow(dX)`.
#' @param Sigma Covariance matrix of size `ncol(dX) x ncol(dX)`.
#' @return A residual matrix of the same size as `dX`.
#' @details The residuals are calculated as
#' \deqn{
#' Z = T^{-1/2} (\Delta X - \mu \Delta T) \Sigma^{-1/2}
#' }{
#' Z = T^{-1/2} (\Delta X - \mu \Delta T) \Sigma^{-1/2}
#' }
#' where \eqn{T} is the Toeplitz matrix whose first column is `acf`.
#'
#' The "square roots" correspond to the Cholesky decomposition for the Toeplitz row-dependence matrix, and the eigen decomposition for the column-dependence matrix.
#'
#' @example examples/fit_setup.R
#' @example examples/csi_resid.R
#'
#' @export
csi_resid <- function(dX, dT, mu, acf, Sigma) {
  # time decorrelation
  Z <- t(t(dX) - mu*dT)
  Z <- cholXZ(X = Z, acf = acf)
  # spatial decorrelation
  ed <- eigen(Sigma)
  C <- t(t(ed$vec) * sqrt(1/ed$val))
  Z %*% C
}