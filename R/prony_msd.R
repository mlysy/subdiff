#' MSD of the Prony-GLE model.
#'
#' @param t Vector of timepoints at which to calculate the MSD.
#' @template args-lambda
#' @template args-nu
#' @param ... Additional arguments to [prony_coeff()].
#' @return Vector of mean square displacements at the timepoints in `t`.
#'
#' @template details-prony_gle
#' @note Since temperature is not provided, the result is only proportional to the desired MSD, such that
#' ```
#' MSD_X(t) = k_B T * prony_msd(t, lambda, nu)
#' ```
#' @export
prony_msd <- function(t, lambda, nu = 1, ...) {
  N <- length(t) # length of output
  msd <- rep(0, N)
  # prony coefficients
  rC <- prony_coeff(lambda, nu, ...)
  r <- rC$r
  C <- rC$C
  K <- length(C) # number of modes
  # msd calculations
  if(K > 1) {
    B <- C[2:K]^2/r
    B0 <- sum(B)
    for(ii in 1:(K-1)) {
      msd <- msd + B[ii] * exp(-r[ii]*abs(t))
    }
  } else B0 <- 0
  C[1]^2 * t + B0 - msd
}
