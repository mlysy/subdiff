#' Fit the farma model.
#'
#' Fit a farma(p,q) model to a multi-dimensional CSI process (See **Details**).
#'
#' @template args-dX
#' @template args-dt
#' @param order A specification of the farma model: the two integer components (p, q) are the AR order and the MA order.
#' @template args-Tz
#' @template args-vcov
#' @return A vector of estimated parameters on transformed scale (See [farma_model()]). If `vcov`, a list with components:
#' \describe{
#' \item{coef}{A vector of estimated parameters on transformed scale.}
#' \item{vcov}{A matrix of estimated covariance of parameters on transformed scale.}
#' }
#'
#' @details The farma(p,q) model is of following form:
#' \deqn{
#' Y_n = \sum_{i=1}^p \phi_i Y_{n-i} + \sum_{j=0}^q \rho_j X_{n-j}
#' }{
#' Y[n] = \phi_1 Y_(n-1) + ... + \phi_p Y_(n-p) + \rho_0 X_(n) + ... + \rho_q X_(n-q)
#' }
#' where \eqn{X_n} is an `q = 1,2` dimensional fBM model (See [fbm_fit()]).
#'
#' Optimization is done by `stats::optimize`. It works better when parameters are re-parametrized into unrestricted form (See [farma_model()]).
#'
#' @example examples/fit_setup.R
#' @example examples/farma_fit.R
#'
#' @export
farma_fit <- function(dX, dt, order, Tz, vcov = TRUE) {
  model <- farma_model(order[1], order[2])
  ans <- csi_fit(model, dX, dt, Tz, vcov)
  ans
}
