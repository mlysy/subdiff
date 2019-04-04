#' Calculate the autocorrelation function for the fBM model.
#' 
#' Compute the autocorrelation of fBM increments of equally-spaced observations (see \strong{Details}).
#' 
#' @param alpha Power law exponent.  A scalar between 0 and 2.
#' @template args-dT
#' @template args-N
#' 
#' @template ret-acf
#' 
#' @details The autocorrelation of fBM increments at lag \eqn{h} is given by
#' \deqn{
#' \textrm{acf}(h) = \Delta t^\alpha/2 \times \left(|n+1|^\alpha + |n-1|^\alpha - 2|n|^\alpha\right).
#' }{
#' acf(n) = \Delta t^\alpha/2 * (|n+1|^\alpha + |n-1|^\alpha - 2|n|^\alpha). 
#' }
#'
#' @example examples/acf_setup.R
#' @example examples/fbm_acf.R
#' @export
fbm_acf <- function(alpha, dT, N) {
  if(N == 1) {
    acf <- dT^alpha
  } else {
    acf <- (dT*(0:N))^alpha
    acf <- .5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2*acf[1:N])
  }
  acf
}
