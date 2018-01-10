#' Model-based residual calculations.
#'
#' Extracts the white-noise residuals from the Gaussian models provided by the package.
#'
#' @param theta Vector of parameter values on the transformed scale returned by the corresponding \code{model_fit} procedure (see Details).
#' @param dX One or two-column matrix of increments.
#' @param dT Interobservation time.
#' @param type,ds Additional arguments specifying for the downsampling estimator (see \code{\link{fds_fit}}).
#' @return A one or two-column matrix of white noise residuals (see Details).
#' @details Each of the models provided by the package generates a normal data matrix \eqn{\Delta X} via the linear transformation
#' \deqn{
#' \Delta X = A Z B + C,
#' }
#' where the variance matrices \eqn{V = AA'} and \eqn{\Sigma = B'B} and mean matrix \eqn{C} are defined by the model parameters, and \eqn{Z} is a matrix of iid \eqn{N(0,1)} white-noise residuals.  The specific square-roots of the variance matrices are the Cholesky decomposition of \eqn{V} (time-domain variance), and the eigendecomposition of \eqn{\Sigma} (spatial-domain variance).
#' @name subdiff-resid

#' @rdname subdiff-resid
#' @export
fbm_resid <- function(theta, dX, dT) {
  q <- ncol(dX) # problem dimensions
  N <- nrow(dX)
  nq <- if(q == 1) 1 else 3
  alpha <- itrans_alpha(theta[1]) # parameters
  mu <- theta[1+1:q]
  Sigma <- itrans_Sigma(theta[q+1+1:nq])
  lsc_resid(dX, dT, mu, fbm_acf(alpha, dT, N), Sigma)
}

#' @rdname subdiff-resid
#' @export
fma_resid <- function(theta, dX, dT) {
  q <- ncol(dX) # problem dimensions
  N <- nrow(dX)
  nq <- if(q == 1) 1 else 3
  alpha <- itrans_alpha(theta[1]) # parameters
  rho <- itrans_rho(theta[2])
  mu <- theta[2+1:q]
  Sigma <- itrans_Sigma(theta[q+2+1:nq])
  dY <- ma_resid(dX, rho)
  lsc_resid(dY, dT, mu, fbm_acf(alpha, dT, N), Sigma)
}

#' @rdname subdiff-resid
#' @export
far_resid <- function(theta, dX, dT) {
  q <- ncol(dX) # problem dimensions
  N <- nrow(dX)
  nq <- if(q == 1) 1 else 3
  alpha <- itrans_alpha(theta[1]) # parameters
  rho <- itrans_rho(theta[2])
  mu <- theta[2+1:q]
  Sigma <- itrans_Sigma(theta[q+2+1:nq])
  dY <- ar_resid(dX, rho)
  lsc_resid(dY, dT, mu, fbm_acf(alpha, dT, N), Sigma)
}

#' @rdname subdiff-resid
#' @export
fds_resid <- function(theta, dX, dT, type, ds) {
  Xt <- apply(rbind(0, dX), 2, cumsum) # recover Xt
  q <- ncol(Xt) # problem dimensions
  N <- nrow(Xt)
  N_ds <- floor(N/ds) - 1
  nq <- if(q == 1) 1 else 3
  alpha <- itrans_alpha(theta[1]) # parameters
  mu <- theta[1+1:q]
  Sigma <- itrans_Sigma(theta[1+q+1:nq])
  if(type == "naive") {
    Xt_ds <- downSample(Xt, ds)
    dX_ds <- apply(Xt_ds, 2, diff)
    ans <- lsc_resid(dX_ds, ds*dT, mu, fbm_acf(alpha, ds*dT, N_ds), Sigma)
  } else {
    dX <- apply(Xt, 2, diff)
    ans <- lsc_resid(dX, dT, mu, fbm_acf(alpha, dT, N-1), Sigma)
  }
}
