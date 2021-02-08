#' Fit the fBM model.
#'
#' Fit a fBM model to a multi-dimensional CSI process.
#'
#' @template args-dX
#' @template args-dt
#' @template args-drift_preset
#' @template args-vcov
#' @template args-ad_only
#' @return A vector of estimated parameters on transformed scale (See [fbm_model()]). If `vcov`, a list with components:
#' \describe{
#' \item{coef}{A vector of estimated parameters on transformed scale.}
#' \item{vcov}{A matrix of estimated covariance of parameters on transformed scale.}
#' }
#'
#' @details The fBM model is of the form:
#' \deqn{
#' X_t = \mu t + \Sigma^{1/2} Z_t
#' }{
#' X[t] = \mu t + \Sigma^{1/2} Z[t],
#' }
#' where \eqn{Z[t]} consists of `q = 1,2` iid fBM processes with \eqn{MSD_Z(t) = t^\alpha}.
#'
#' Optimization is done by [stats::optimize()]. It works better when parameters are re-parametrized into unrestricted form (See [fbm_model()]).
#'
#' @example examples/fit_setup.R
#' @example examples/fbm_fit.R
#'
#' @export
fbm_fit <- function(dX, dt, drift = c("linear", "none", "quadratic"),
                    vcov = TRUE, ad_only = TRUE) {
  ## model <- fbm_model()
  ## ans <- csi_fit(model, dX, dt, Tz, vcov)
  ## ans
  drift <- match.arg(drift)
  model <- fbm_model$new(dX = dX, dt = dt, drift = drift)
  out <- model$fit(psi0 = c(-5, 5), vcov = vcov)
  if(ad_only) out <- to_ad(out, model = model, vcov = vcov)
  out
}

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
