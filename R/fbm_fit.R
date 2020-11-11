#' Fit the fBM model.
#'
#' Fit a fBM model to a multi-dimensional CSI process.
#'
#' @template args-dX
#' @template args-dt
#' @template args-drift_preset
#' @template args-vcov
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
#' Optimization is done by [optim()]. It works better when parameters are re-parametrized into unrestricted form (See [fbm_model()]).
#'
#' @example examples/fit_setup.R
#' @example examples/fbm_fit.R
#'
#' @export
fbm_fit <- function(dX, dt, drift = c("linear", "none", "quadratic"),
                    vcov = TRUE) {
  ## model <- fbm_model()
  ## ans <- csi_fit(model, dX, dt, Tz, vcov)
  ## ans
  mobj <- fbm_model$new(dX = dX, dt = dt, drift = drift)
  mobj$fit(phi0 = c(.001, 1.999), vcov = vcov)
}
