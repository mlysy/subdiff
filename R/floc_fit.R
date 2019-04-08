#' Fit the floc model.
#' 
#' Fit the Savin & Doyle's localization error model with fBM process (See \strong{Details}).
#'
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @return A vector of estimated parameters on transformed scale (See \code{\link{floc_model}}). If \code{var_calc}, a list with components:
#' \describe{
#' \item{coef}{A vector of estimated parameters on transformed scale.}
#' \item{vcov}{A matrix of estimated covariance of parameters on transformed scale.}
#' }
#'
#' @details The Savin & Doyle's localization error model with fBM process (floc model) has the following form:
#' \deqn{
#' Y_n = 1/\tau \int_0^\tau X_{n \Delta t +s} ds + \sigma e_{n}
#' }{
#' Y[n] = 1/\tau integral_0^\tau X_{n  \Delta t +s} ds + \sigma e[n]
#' }
#' where \eqn{X_n} is an \code{q = 1,2} dimensional fBM model (See \code{\link{fbm_fit}}).
#'
#' The integral part measures the average displacement of particle trajectories during \eqn{[n \Delta t, n\Delta t + \tau]} and this term is called dynamic error.
#' \eqn{e_n} is a standard white noise, and this term is called static error.
#' 
#' When put \code{tau} = 0, this model becomes fractional static error model.
#' When put \code{sigma2} = 0, this model becomes fractional dynamic error model.
#' When put \code{(tau, sigma2)} = 0, this model becomes fractional Brownian model.
#' 
#' Optimization is done by \code{\link{optimize}}. It works better when parameters are re-parametrized into unrestricted form (See \code{\link{floc_model}}).
#' 
#' @references Savin, Thierry, and Patrick S. Doyle. "Static and dynamic errors in particle tracking microrheology." Biophysical journal 88.1 (2005): 623-638.
#' 
#' @example examples/fit_setup.R
#' @example examples/floc_fit.R
#'
#' @export
floc_fit <- function(dX, dT, Tz, var_calc = TRUE) {
  model <- floc_model()
  ans <- csi_fit(model, dX, dT, Tz, var_calc)
  ans
}

# #' @param tau Scalar between 0 and 1, indicating the percentage of time for which the camera shutter is open in the dynamic error model.  Estimated if missing. See Details.
# #' @param sigma2 Magnitude of static error.  Estimated if missing. See Details.
