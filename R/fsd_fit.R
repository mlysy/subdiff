#' Fit the fsd model.
#'
#' Fit the Savin & Doyle's localization error model with fBM process (See **Details**).
#'
#' @template args-dX
#' @template args-dt
#' @template args-Tz
#' @template args-vcov
#' @return A vector of estimated parameters on transformed scale (See [fsd_model()]). If `vcov`, a list with components:
#' \describe{
#' \item{coef}{A vector of estimated parameters on transformed scale.}
#' \item{vcov}{A matrix of estimated covariance of parameters on transformed scale.}
#' }
#'
#' @details The Savin & Doyle's localization error model with fBM process (fsd model) has the following form:
#' \deqn{
#' Y_n = 1/\tau \int_0^\tau X_{n \Delta t +s} ds + \sigma e_{n}
#' }{
#' Y[n] = 1/\tau integral_0^\tau X_{n  \Delta t +s} ds + \sigma e[n]
#' }
#' where \eqn{X_n} is an `q = 1,2` dimensional fBM model (See [fbm_fit()]).
#'
#' The integral part measures the average displacement of particle trajectories during \eqn{[n \Delta t, n\Delta t + \tau]} and this term is called dynamic error.
#' \eqn{e_n} is a standard white noise, and this term is called static error.
#'
#' When put `tau` = 0, this model becomes fractional static error model.
#' When put `sigma2` = 0, this model becomes fractional dynamic error model.
#' When put `(tau, sigma2)` = 0, this model becomes fractional Brownian model.
#'
#' Optimization is done by [optimize()]. It works better when parameters are re-parametrized into unrestricted form (See [fsd_model()]).
#'
#' @references Savin, Thierry, and Patrick S. Doyle. "Static and dynamic errors in particle tracking microrheology." Biophysical journal 88.1 (2005): 623-638.
#'
#' @example examples/fit_setup.R
#' @example examples/fsd_fit.R
#'
#' @export
fsd_fit <- function(dX, dt, Tz, vcov = TRUE) {
  model <- fsd_model()
  ans <- csi_fit(model, dX, dt, Tz, vcov)
  ans
}

# #' @param tau Scalar between 0 and 1, indicating the percentage of time for which the camera shutter is open in the dynamic error model.  Estimated if missing. See Details.
# #' @param sigma2 Magnitude of static error.  Estimated if missing. See Details.
