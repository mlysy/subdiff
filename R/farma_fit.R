#' Fit the ARMA(p,q) filter with fBM residuals.
#'
#' @name farma_fit
#' @template args-dX
#' @template args-dT
#' @param order A specification of the ARMA filter: the two integer components (p, q) are the AR order and the MA order.
#' @template args-Tz
#' @template args-var_calc
#' @template ret-cov_vcov
#'
#' @details The ARMA(p,q) filter is of the form
#' \deqn{
#' Y_t = \sum_{i=1}^p \phi_i Y_{t-i} + \sum_{j=0}^q \rho_j X_{t-j}, \rho_0 = 1-\sum_{i=1}^p \phi_i-\sum_{j=1}^q \rho_j,
#' }
#' where residuals \eqn{X_t} follows the fBM model:
#' \deqn{
#' X_t = \mu t + \Sigma^{1/2} Z_t,
#' }
#' where \eqn{Z_t} consists of \code{q = 1,2} iid fBM processes.
#'
#' @example examples/parameters.R
#' @examples
#' fits <- farma_fit(dX, dT, order = c(1,1), Tz, var_calc = TRUE) # fitting ARMA(1,1) filter
#'
#' @export
farma_fit <- function(dX, dT, order, Tz, var_calc = TRUE) {
  model <- farma_model(order[1], order[2])
  ans <- csi_fit(model, dX, dT, Tz, var_calc)
  ans
}
