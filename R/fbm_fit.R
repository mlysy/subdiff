#' Fit the fBM model.
#'
#' @template args-Xt
#' @template args-dt
#' @template args-drift_preset
#' @template args-vcov
#' @template args-ad_only
#'
#' @template ret-fit
#'
#' @seealso [fbm_model], the class definition for the fBM model.
#'
#' @example examples/fbm_sim.R
#' @example examples/fbm_fit.R
#'
#' @export
fbm_fit <- function(Xt, dt, drift = c("linear", "none", "quadratic"),
                    vcov = TRUE, ad_only = TRUE) {
  ## model <- fbm_model()
  ## ans <- csi_fit(model, dX, dt, Tz, vcov)
  ## ans
  drift <- match.arg(drift)
  model <- fbm_model$new(Xt = Xt, dt = dt, drift = drift)
  out <- model$fit(psi0 = c(-5, 5), vcov = vcov)
  if(ad_only) out <- to_ad(out, model = model, vcov = vcov)
  out
}

## ' @details The fBM model is of the form:
## ' \deqn{
## ' X_t = \mu t + \Sigma^{1/2} Z_t
## ' }{
## ' X[t] = \mu t + \Sigma^{1/2} Z[t],
## ' }
## ' where \eqn{Z[t]} consists of `q = 1,2` iid fBM processes with \eqn{MSD_Z(t) = t^\alpha}.
## '
## ' Optimization is done by [stats::optimize()]. It works better when parameters are re-parametrized into unrestricted form (See [fbm_model()]).
## '


#--- helper functions ----------------------------------------------------------

#' Convert an estimate of `omega` to an estimate of `ad`.
#'
#' @noRd
to_ad <- function(omega, model, vcov) {
  if(!vcov) {
    ad_hat <- model$get_subdiff(omega)
    return(ad_hat)
  } else {
    ad_hat <- model$get_subdiff(omega$coef)
    # delta method
    Jac <- numDeriv::jacobian(func = model$get_subdiff, x = omega$coef)
    ad_vcov <- tcrossprod(Jac %*% omega$vcov, Jac)
    colnames(ad_vcov) <- rownames(ad_vcov) <- names(ad_hat)
    return(list(coef = ad_hat, vcov = ad_vcov))
  }
}
