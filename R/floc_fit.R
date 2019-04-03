#' Fit the fBM model with Savin & Doyle's Localization Error.
#'
#' @name floc_fit
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @template ret-cov_vcov 
#' 
#' @details The fBM + localization error (fLOC) model has the form
#' \deqn{
#' Y_n = 1/\tau \int_0^\tau X_{n+s} ds + \sigma e_{n},
#' },
#' where \eqn{X_n} is a 1D or 2D fBM process
#' \deqn{
#' X_t = \mu t + \Sigma^{1/2} Z_t,
#' }. 
#' where \eqn{Z_t} consists of \code{q = 1,2} iid fBM processes.
#' When put \code{tau} = 0, this model becomes fractional localization error model.
#' When put \code{sigma2} = 0, this model becomes fractional dynamic error model.
#' When put \code{(tau, sigma2)} = 0, this model becomes fractional Brownian model.
#' 
#' @example example/parameters.R
#' @examples 
#' fits <- floc_fit(dX, dT, Tz, var_calc = TRUE)
#' 
#' @export
floc_fit <- function(dX, dT, Tz, var_calc = TRUE) {
  model <- floc_model()
  ans <- csi_fit(model, dX, dT, Tz, var_calc)
  ans
}

# #' @param tau Scalar between 0 and 1, indicating the percentage of time for which the camera shutter is open in the dynamic error model.  Estimated if missing. See Details.
# #' @param sigma2 Magnitude of static error.  Estimated if missing. See Details.
