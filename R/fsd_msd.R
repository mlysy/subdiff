#' Calculate the mean square displacement of the fSD model.
#'
#' Compute the MSD of the Savin & Doyle (2005) localization error model with fBM increments (see [fsd_acf()]).
#'
#' @template args-t
#' @template args-dt
#' @param alpha Subdiffusion exponent of the underlying fBM process. A scalar between 0 and 2.
#' @param tau The ratio between camera shutter open time and the interobservation time `dt`. A scalar between 0 and 1.
#' @param sigma2 The magnitude of the static error. A positive scalar.
#'
#' @template ret-msd
#'
#' @template details-fsd_model
#'
#' @template ref-savin_doyle
#'
#' @export
fsd_msd <- function(t, dt, alpha, tau, sigma2) {
  if(tau == 0) {
    msd <- fbm_msd(t, alpha)
  } else {
    tau2 <- tau * dt # convert to absolute units
    msd <- abs(t + tau2)^(alpha+2) + abs(t - tau2)^(alpha+2) - 2 * abs(t)^(alpha+2)
    msd <- msd - 2 * abs(tau2)^(alpha+2)
    msd <- msd / (tau2^2 * (alpha+1) * (alpha+2))
  }
  msd <- msd + 2 * sigma2
  msd
}
