#' Autocorrelation of the Prony-GLE increment process.
#'
#' @template args-lambda
#' @template args-nu
#' @template args-dt
#' @template args-N
#' @param ... Additional arguments to [prony_coeff()].
#'
#' @template ret-acf
#' @template details-prony_gle
#' @note Since temperature is not provided, the result is only proportional to the desired ACF, such that
#' ```
#' ACF_dX(t) = k_B T * prony_acf(lambda, nu, N, dt)
#' ```
#' @export
prony_acf <- function(lambda, nu = 1, dt, N, ...) {
  acf <- rep(0, N)
  # prony coefficients
  rC <- prony_coeff(lambda, nu, ...)
  r <- rC$r
  C <- rC$C
  K <- length(C) # number of modes
  # acf calculations
  if(K > 1) {
    for(ii in 1:(K-1)) {
      tmp2 <- exp(-r[ii]*(0:N)*dt)
      tmp3 <- 2*tmp2[1:N] - tmp2[1:N+1] - c(tmp2[2], tmp2[1:(N-1)])
      acf <- acf + C[ii+1]^2/(2*r[ii]) * tmp3
    }
  }
  acf[1] <- acf[1] + C[1]^2*dt
  acf
}
