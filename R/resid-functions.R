#' Model-based residual calculations.
#'
#' Extracts the white-noise residuals from the Gaussian models provided by the package.
#'
#' @param theta Vector of parameter values on the transformed scale returned by the corresponding `model_fit` procedure (see Details).
#' @param dX One or two-column matrix of increments.
#' @param dt Interobservation time.
#' @return A one or two-column matrix of white noise residuals (see Details).
#' @details Each of the models provided by the package generates a normal data matrix \eqn{\Delta X} via the linear transformation
#' \deqn{
#' \Delta X = A Z B + C,
#' }
#' where the variance matrices \eqn{V = AA'} and \eqn{\Sigma = B'B} and mean matrix \eqn{C} are defined by the model parameters, and \eqn{Z} is a matrix of iid \eqn{N(0,1)} white-noise residuals.  The specific square-roots of the variance matrices are the Cholesky decomposition of \eqn{V} (time-domain variance), and the eigendecomposition of \eqn{\Sigma} (spatial-domain variance).
#' @name subdiff-resid

#' @rdname subdiff-resid
#' @export
fbm_resid <- function(theta, dX, dt) {
  qq <- ncol(dX) # problem dimensions
  N <- nrow(dX)
  nq <- if(qq == 1) 1 else 3
  alpha <- itrans_alpha(theta[1]) # parameters
  mu <- theta[1+1:qq]
  Sigma <- itrans_Sigma(theta[qq+1+1:nq])
  csi_resid(dX, dt, mu, fbm_acf(alpha, dt, N), Sigma)
}

#' @rdname subdiff-resid
#' @param order Integer vector of length 2 specifying the number of AR and MA terms respectively.
#' @export
farma_resid <- function(theta, dX, dt, order) {
  qq <- ncol(dX) # problem dimensions
  N <- nrow(dX)
  nq <- if(qq == 1) 1 else 3
  p <- order[1] # ARMA order
  q <- order[2]
  alpha <- itrans_alpha(theta[1]) # parameters
  if(p) {
    phi <- theta[1+1:p]
  } else {
    phi <- 0
  }
  rho <- theta[1+p+1:q]
  mu <- theta[1+p+q+1:qq]
  Sigma <- itrans_Sigma(theta[qq+1+p+q+1:nq])
  csi_resid(dX, dt, mu, farma_acf(alpha, phi, rho, dt, N), Sigma)
}

#' @rdname subdiff-resid
#' @export
floc_resid <- function(theta, dX, dt) {
  qq <- ncol(dX) # problem dimensions
  N <- nrow(dX)
  nq <- if(qq == 1) 1 else 3
  alpha <- itrans_alpha(theta[1]) # parameters
  tau <- itrans_tau(theta[2])
  sigma2 <- exp(2*theta[3])
  mu <- theta[3+1:qq]
  Sigma <- itrans_Sigma(theta[qq+3+1:nq])
  acf1 <- floc_acf(alpha, tau, sigma2, dt, N)
  res <- csi_resid(dX, dt, mu, acf1, Sigma)
  res
}

# fdl_resid <- function(theta, dX, dt, type = c("fdl", "fdy", "flo")) {
#   qq <- ncol(dX) # problem dimensions
#   N <- nrow(dX)
#   nq <- if(qq == 1) 1 else 3
#   type <- match.arg(type)
#   if(type == "fdl") {
#     alpha <- itrans_alpha(theta[1]) # parameters
#     tau <- itrans_tau(theta[2])
#     sigma2 <- exp(2*theta[3])
#     mu <- theta[3+1:qq]
#     Sigma <- itrans_Sigma(theta[qq+3+1:nq])
#     acf1 <- fdyn_acf(alpha, tau, dt, N) + sigma2 * c(2, 1, rep(0, N-2))
#     res <- csi_resid(dX, dt, mu, acf1, Sigma)
#   } else if(type == "fdy") {
#     alpha <- itrans_alpha(theta[1]) # parameters
#     tau <- itrans_tau(theta[2])
#     mu <- theta[2+1:qq]
#     Sigma <- itrans_Sigma(theta[qq+2+1:nq])
#     acf1 <- fdyn_acf(alpha, tau, dt, N)
#     res <- csi_resid(dX, dt, mu, acf1, Sigma)
#   } else if(type == "flo") {
#     alpha <- itrans_alpha(theta[1]) # parameters
#     sigma2 <- theta[2]^2
#     mu <- theta[2+1:qq]
#     Sigma <- itrans_Sigma(theta[qq+2+1:nq])
#     acf1 <- fbm_acf(alpha, dt, N) + sigma2 * c(2, 1, rep(0, N-2))
#     res <- csi_resid(dX, dt, mu, acf1, Sigma)
#   }
#   res
# }
