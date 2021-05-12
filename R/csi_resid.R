#' Calculate the residuals of a CSI process.
#'
#' @template args-dX
#' @param drift Matrix of the same size as `dX` of increment drift in each coordinate per timepoint (see 'Details').
#' @param acf Autocorrelation vector of length `nrow(dX)` (see 'Details').
#' @param Sigma Covariance matrix of size `ncol(dX) x ncol(dX)`.
#' @return A residual matrix of the same size as `dX`.
#' @details The residuals are calculated as
#' ```
#' Z = toeplitz(acf)^{-1/2} %*% (dX - drift) %*% Sigma^{-1/2},
#' ```
#' where the "square roots" correspond to the Cholesky decomposition for the row-dependence matrix `toeplitz(acf)`, and the eigen decomposition for the column-dependence matrix `Sigma`.
#'
#' @example examples/fbm_sim.R
#' @example examples/csi_resid.R
#'
#' @export
csi_resid <- function(dX, drift, acf, Sigma) {
  # time decorrelation
  ## Z <- t(t(dX) - mu*dt)
  Z <- cholXZ(X = dX - drift, acf = acf)
  # spatial decorrelation
  ed <- eigen(Sigma)
  C <- t(t(ed$vectors) * sqrt(1/ed$values))
  Z %*% C
}
