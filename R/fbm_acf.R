#' Calculate the ACF for the fBM model.
#'
#' Compute the autocovariance of fBM increments of equally-spaced observations (see **Details**).
#'
#' @param alpha Power law exponent. A scalar between 0 and 2.
#' @template args-dt
#' @template args-N
#'
#' @template ret-acf
#'
#' @details The autocovariance of fBM increments at lag \eqn{h} is given by
#' \deqn{
#' \textrm{acf}(h) = \Delta t^\alpha/2 \times \left(|h+1|^\alpha + |h-1|^\alpha - 2|h|^\alpha\right).
#' }{
#' acf(h) = \Delta t^\alpha/2 * (|h+1|^\alpha + |h-1|^\alpha - 2|h|^\alpha).
#' }
#'
#' This function returns the autocovariance function of increments of fBM process.
#'
#' @example examples/acf_setup.R
#' @example examples/fbm_acf.R
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
