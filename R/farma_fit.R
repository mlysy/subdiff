#' Fit the farma model.
#' 
#' Fit a farma(p,q) model to a multi-dimensional CSI process (See \strong{Details}).
#'
#' @template args-dX
#' @template args-dT
#' @param order A specification of the farma model: the two integer components (p, q) are the AR order and the MA order.
#' @template args-Tz
#' @template args-var_calc
#' @return A vector of estimated parameters on transformed scale (See \code{\link{farma_model}}). If \code{var_calc}, a list with components:
#' \describe{
#' \item{coef}{A vector of estimated parameters on transformed scale.}
#' \item{vcov}{A matrix of estimated covariance of parameters on transformed scale.}
#' }
#'
#' @details The farma(p,q) model is of following form:
#' \deqn{
#' Y_n = \sum_{i=1}^p \phi_i Y_{n-i} + \sum_{j=0}^q \rho_j X_{n-j}
#' }{
#' Y[n] = \phi_1 Y[n-1] + ... + \phi_p Y[n-p] + \rho_0 X[n] + ... + \rho_q X[n-q]
#' }
#' where \eqn{X_n} is an \code{q = 1,2} dimensional fBM model (See \code{\link{fbm_fit}}).
#' 
#' Optimization is done by \code{\link{optimize}}. It works better when parameters are re-parametrized into unrestricted form (See \code{\link{farma_model}}).
#'
#' @example examples/fit_setup.R
#' @example examples/farma_fit.R
#'
#' @export
farma_fit <- function(dX, dT, order, Tz, var_calc = TRUE) {
  model <- farma_model(order[1], order[2])
  ans <- csi_fit(model, dX, dT, Tz, var_calc)
  ans
}
