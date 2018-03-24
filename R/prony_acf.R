#' @title Prony Acf.
#' @description GLE with 0 mass and 0 external potential.
#' @details
#' process \code{x} is defined as
#' \code{vsigma * int_{-infty}^t gamma(t-s)  x_s d s = k_B T vsigma * F_t}
#' where the force noise F_t is a sum of independent OU processes, ie,
#' cov(F_s, F_{t+s}) = gamma(t) = exp(-alpha1 t) + ... + exp(-alphaK t).
#' the resulting process is
#' Y_t = C1 X_1t + ... + CK X_Kt,
#' where X_1t = B_1t and X_it, i > 1 are OU processes,
#' d X_it = -ri X_it dt + d B_it,
#' and all B_it are independent.
#' note that cov(X_i0, X_it) = (2*ri)^{-1} * exp(-ri * t).
#' @param lambda coefficients of the sum of OU processes.
#' @param N number of samples.
#' @param dt interobservation time.
#' @param ... Additional arguments to \code{\link{prony_coeff}} for obtaining the BM-OU coefficients.
#' @return Vector of autocorrelations.
#' @export
prony_acf <- function(lambda, N, dt, ...) {
  acf <- rep(0, N)
  # prony coefficients
  rC <- prony_coeff(lambda, ...)
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
