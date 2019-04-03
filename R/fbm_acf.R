#' ACF of the increment of fBM process (fractional Gaussian noise).
#' 
#' @template args-alpha
#' @template args-dT
#' @template args-N
#' @template ret-acf
#' 
#' @details The expression of the ACF of fractional Gaussian noise is
#' \deqn{
#' \frac{1}{2}(|n-1|^\alpha + |n + 1|^\alpha - 2|n|^\alpha) \Delta t^\alpha
#' }
#' 
#' @examples 
#' acf1 <- fbm_acf(alpha = 0.8, dT = 1/60, N = 1800)
#' 
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
