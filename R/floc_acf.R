#' ACF of fBM + Savin & Doyle's Localization Error.
#'
#' @template args-alpha
#' @param tau Real number ranging between 0 and 1, indicating the percentage of time for which the camera shutter is open in the dynamic error model. See Details.
#' @param sigma2 Magnitude of static error. See Details.
#' @template args-dT
#' @template args-N
#' @template ret-acf
#' 
#' @details this function returns the ACF of \eqn{\Delta Y_n}, where \eqn{Y_n} follows fLOC model:
#' \deqn{
#' Y_n = 1/\tau \int_0^\tau X_{n+s} ds + \sigma e_{n},
#' },
#' where \eqn{X_n} is an fBM process and \eqn{e_n} is a white noise with variance \code{sigma2}.
#' 
#' @examples 
#' floc_acf(alpha = 0.8, tau = 1/60, sigma2 = 0.01, dT = 1/60, N = 1800)
#' 
#' @export
floc_acf <- function(alpha, tau, sigma2, dT, N) {
  if(!tau) {
    acf1 <- fbm_acf(alpha, dT, N)
  } else {
    vec <- .g_func(alpha, tau, dT, N)
    if(N == 1) {
      acf1 <- 2 * vec[2] - 2 * vec[1]
    } else {
      acf1 <- vec[1:N + 1] + c(vec[2], vec[1:(N - 1)]) - 2 * vec[1:N]  
    }
  }
  acf1[1:2] <- acf1[1:2] + sigma2*c(2,-1)
  acf1
}

.g_func <- function(alpha, tau, dT, N) {
  vec <- (0:N) / tau
  alpha2 <- alpha + 2
  ans <- (vec + 1)^alpha2 + abs(vec - 1)^alpha2 - 2 * vec^alpha2
  ans * (tau * dT)^alpha / (alpha+1) / alpha2 / 2
}

