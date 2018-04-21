#' Residuals of a location-scale stationary increments model.
#'
#' @param dX Matrix of increments where each column is a time series.
#' @param dT Interobservation time.
#' @param mu Vector of linear drift in each direction.
#' @param acf Autocorrelation, a vector of length \code{nrow(dX)}.
#' @param Sigma Covariance matrix of size \code{ncol(dX) x ncol(dX)}.
#' @return A residual matrix of the same size as \code{dX}.
#' @details The residuals are calculated as
#' \deqn{
#' Z = Toeplitz(acf)^{-1/2} (dX - mu dT) Sigma^{-1/2}
#' }
#' where the "square roots" correspond to the Cholesky decomposition for the Toeplitz row-dependence matrix, and the eigen decomposition for the column-dependence matrix.
#' @export
lsc_resid <- function(dX, dT, mu, acf, Sigma) {
  # time decorrelation
  Z <- t(t(dX) - mu*dT)
  Z <- cholXZ(X = Z, acf = acf)
  # spatial decorrelation
  ed <- eigen(Sigma)
  C <- t(t(ed$vec) * sqrt(1/ed$val))
  Z %*% C
}
