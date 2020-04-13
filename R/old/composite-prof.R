#' @title Profile Composite Likelihood
#' @param Y Original time series following Matrix Normal MN(X*Beta, V, Sigma), downsampled case 
#' would be Y_ds ~ MN(X_ds*Beta, V_ds, Sigma)
#' @param X Linear drift of time series. If X is of length 1, X_ds = rep(X, n)
#' @param acf ACF of columnwise-variance matrix `V_ds`, either vector or Toeplitz-object
#' @param ds Downsampling rate, starting at 2
#' @param noSigma logic, assume Sigma is identity matrix if true.
#' @return Composite likelihood likelihood using profile sufficient statistics
#' @export
composite.prof <- function(suff, Y, X, acf, ds, noSigma = FALSE) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- composite.suff(Y = Y, X = X, acf = acf, ds = ds)
  }
  n.ds <- suff$n.ds
  S <- suff$S
  ldV <- suff$ldV
  q <- nrow(S)
  if(!noSigma) {
    ll <- n.ds*determinant(S, logarithm = TRUE)$mod[1] + q*ldV
  } else {
    ll <- sum(diag(S)) + q*ldV
  }
  -.5 * ds * ll
}
