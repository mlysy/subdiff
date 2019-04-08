#' Fit the fBM model.
#'
#' Fit a fBM model to a multi-dimensional CSI process.
#'
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @return A vector of estimated parameters on transformed scale (See \code{\link{fbm_model}}). If \code{var_calc}, a list with components:
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
#' where \eqn{Z[t]} consists of \code{q = 1,2} iid fBM processes with \eqn{MSD_Z(t) = t^\alpha}.
#' 
#' Optimization is done by \code{\link{optim}}. It works better when parameters are re-parametrized into unrestricted form (See \code{\link{fbm_model}}).
#' 
#' @example examples/fit_setup.R
#' @example examples/fbm_fit.R
#'
#' @export
fbm_fit <- function(dX, dT, Tz, var_calc = TRUE) {
  model <- fbm_model()
  ans <- csi_fit(model, dX, dT, Tz, var_calc)
  ans
}
