#' Fit the fBM model.
#' 
#' @name fbm_fit
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @template ret-cov_vcov
#' 
#' @details The fBM model is of the form
#' \deqn{
#' X_t = \mu t + \Sigma^{1/2} Z_t,
#' }
#' where \eqn{Z_t} consists of \code{q = 1,2} iid fBM processes with \eqn{MSD_Z(t) = t^\alpha}.
#' @example example/parameters.R
#' 
#' @export
fbm_fit <- function(dX, dT, Tz, var_calc = TRUE) {
  model <- fbm_model()
  ans <- csi_fit(model, dX, dT, Tz, var_calc)
  ans
}
