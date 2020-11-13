#' Calculate the ACF for the fsd model.
#'
#' Compute the autocovariance for the increments of Savin & Doyle's localization model with fBM process (see **Details**).
#'
#' @param alpha Power law exponent of fBM process. A scalar between 0 and 2.
#' @param tau The ratio between camera shutter open time and interobservation time `dt` (see **Details**). A scalar between 0 and 1.
#' @param sigma2 The Magnitude of static errors (see **Details**). A positive scalar.
#' @template args-dt
#' @template args-N
#'
#' @template ret-acf
#'
#' @details The Savin & Doyle's localization error model with fBM process (fsd model) is of following form:
#' \deqn{
#' Y_n = 1/\tau \int_0^\tau X_{n  \Delta t+s} ds + \sigma e_{n}
#' }{
#' Y[n] = 1/\tau integral_0^\tau X_{n  \Delta t+s} ds + \sigma e[n]
#' }
#' where \eqn{X_n} is an fBM process.
#'
#' The integral part measures the average displacement of particle trajectories during \eqn{[n \Delta t, n\Delta t + \tau]} and this term is called dynamic error.
#' \eqn{e[n]} is a standard white noise, and this term is called static error.
#'
#' This function returns the autocovariance function of increments of \eqn{Y[n]}.
#'
#' @references Savin, T., and Doyle, P.S. "Static and dynamic errors in particle tracking microrheology." Biophysical Journal 88.1 (2005): 623-638.
#'
#' @example examples/acf_setup.R
#' @example examples/fsd_acf.R
#'
#' @export
fsd_acf <- function(alpha, tau, sigma2, dt, N) {
  if(!tau) {
    acf1 <- fbm_acf(alpha, dt, N)
  } else {
    vec <- .g_func(alpha, tau, dt, N)
    if(N == 1) {
      acf1 <- 2 * vec[2] - 2 * vec[1]
    } else {
      acf1 <- vec[1:N + 1] + c(vec[2], vec[1:(N - 1)]) - 2 * vec[1:N]
    }
  }
  acf1[1:2] <- acf1[1:2] + sigma2*c(2,-1)
  acf1
}

.g_func <- function(alpha, tau, dt, N) {
  vec <- (0:N) / tau
  alpha2 <- alpha + 2
  ans <- (vec + 1)^alpha2 + abs(vec - 1)^alpha2 - 2 * vec^alpha2
  ans * (tau * dt)^alpha / (alpha+1) / alpha2 / 2
}

