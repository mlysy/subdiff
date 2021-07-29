#' ACF2 of fBM + Dynamic Error
#'
#' @param alpha Subdiffusion exponent
#' @param tau Ratio of time-window over `dT`, range from 0 to 1.
#' @param dT Vector of time points \eqn{ \{\Delta t, 2\Delta t, ..., N\Delta t\} }
#' @param N length
#' @details 
#' this function returns the MSD of \eqn{Y_t}, the integral of fBM process \eqn{X_t} with subdiffusion 
#' exponent \eqn{\alpha} \deqn{Y_t = \int_{0}^{\tau} X(t-s)ds}. The expression of the MSD is
#' \deqn{\frac{(t+\tau)^\alpha + (t-\tau)^\alpha - 2t^\alpha - 2\tau^\alpha}{(\alpha+1)(\alpha+2)}}
#' @examples 
#' fdyn.msd(alpha = 0.8, tau = 1/600, t = (1:200) * 1/60)
#' @export
fdyn_acf <- function(alpha, tau, dT, N) {
  if(!tau) {
    acf1 <- fbm_acf(alpha, dT, N)
  } else {
    vec <- .g_func2(alpha, tau, dT, N)
    if(N == 1) {
      acf1 <- 2 * vec[2] - 2 * vec[1]
    } else {
      acf1 <- vec[1:N + 1] + c(vec[2], vec[1:(N - 1)]) - 2 * vec[1:N]  
    }
  }
  acf1
}

#' ACF of fBM process with localization errors
#' @export
floc_acf <- function(alpha, tau, sigma2, dT, N) {
  acf1 <- fdyn_acf(alpha, tau, dT, N)
  acf1[1:2] <- acf1[1:2] + sigma2*c(2,-1)
  acf1
}

.g_func2 <- function(alpha, tau, dT, N) {
  vec <- (0:N) / tau
  alpha2 <- alpha + 2
  ans <- (vec + 1)^alpha2 + abs(vec - 1)^alpha2 - 2 * vec^alpha2
  ans * (tau * dT)^alpha / (alpha+1) / alpha2 / 2
}

