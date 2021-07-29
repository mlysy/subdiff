#' Calculate the autocorrelation of the fBM model.
#'
#' Compute the autocorrelation of fractional Brownian motion (fBM) increments at equally-spaced time points (see 'Details').
#'
#' @param alpha Subdiffusion exponent. A scalar between 0 and 2.
#' @template args-dt
#' @template args-N
#'
#' @template ret-acf
#'
#' @details Let `X_t` denote an fBM process and `dX_n = X_{dt * (n+1)} - X_{dt * n}` denote the `n`th increment of `X_t` with interobservation time `dt`.  The autocorrelation of the fBM increment process `dX_n` is given by
#' ```
#' acf_dX(n) = 0.5 * dt^alpha * (|n+1|^alpha + |n-1|^alpha - 2|n|^\alpha).
#' ```
#'
#' @export
fbm_acf <- function(alpha, dt, N) {
  if(N == 1) {
    acf <- dt^alpha
  } else {
    acf <- (dt*(0:N))^alpha
    acf <- .5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2*acf[1:N])
  }
  acf
}
