#' Fit the fARMA(p,q) model.
#'
#' @template args-Xt
#' @template args-dt
#' @param order A specification of the farma model: the two integer components (p, q) are the AR order and the MA order.
#' @template args-drift_preset
#' @template args-vcov
#' @template args-ad_only
#'
#' @template ret-fit
#'
#' @seealso [farma_model], the class definition for the fARMA(p,q) model.
#'
#' @example examples/farma_sim.R
#' @example examples/farma_fit.R
#'
#' @export
farma_fit <- function(Xt, dt, order,
                      drift = c("linear", "none", "quadratic"),
                      vcov = TRUE, ad_only = TRUE) {
  drift <- match.arg(drift)
  ## if((length(order) != 2) || !all(order - as.integer(order) == 0)) {
  ##   stop("order must be a vector of 2 nonnegative integers.")
  ## }
  if(length(order) == 2 && all(order == 0)) {
    return(fbm_fit(Xt = Xt, dt = dt, drift = drift,
                   vcov = vcov, ad_only = ad_only))
  }
  model <- farma_model$new(Xt = Xt, dt = dt, drift = drift,
                           order = order, m = 50)
  out <- model$fit(psi0 = rep(0, length(model$phi_names)),
                   vcov = vcov)
  if(ad_only) out <- to_ad(out, model, vcov)
  out
  ## ans <- csi_fit(model, dX, dt, Tz, vcov)
  ## ans
}

## ' @details The farma(p,q) model is of following form:
## ' \deqn{
## ' Y_n = \sum_{i=1}^p \phi_i Y_{n-i} + \sum_{j=0}^q \rho_j X_{n-j}
## ' }{
## ' Y[n] = \phi_1 Y_(n-1) + ... + \phi_p Y_(n-p) + \rho_0 X_(n) + ... + \rho_q X_(n-q)
## ' }
## ' where \eqn{X_n} is a `d`-dimensional fBM model (See [fbm_fit()]).
## '
## ' Optimization is done by [stats::optim()]. It works better when parameters are re-parametrized into unrestricted form (See [farma_model()]).
## '
