#' MSD of the Prony-GLE model.
#'
#' @details
#' The Prony-GLE model satisfies the integro-differential equation
#' \deqn{
#' F_t - \int_{-\infty}^t \gamma(t-s) \dot X_s ds = 0
#' }
#' \deqn{
#' acf_F(t) = k_B T \gamma(t),
#' }
#' where \eqn{T} is temperature, \eqn{k_B} is Boltzmann's constant, and the memory kernel is a sum of exponentials:
#' \deqn{
#' \gamma(t) = 1/K sum_{n=1}^K exp(-\lambda_n * t).
#' }
#' The solution process is of the form
#' \deqn{
#' X_t = C_0 B_t + \sum_{i=1}^{K-1} C_i W_{it},
#' }
#' where \eqn{d W_{it} = -\rho_i W_{it} dt + d B_{it}} are Ornstein-Uhlenbeck processes all independent of each other and of \eqn{B_t}.
#' @note Since temperature is not provided, the result is only proportional to the desired GLE, such that PUT IN CORRECT EQUATION.
#' @param t times at which to calculate the MSD.
#' @param lambda coefficients of the sum of OU processes.
#' @param ... Additional arguments to \code{prony.coeff} for obtaining the BM-OU coefficients.
#' @return Vector of mean square displacements.
#' @export
prony.msd <- function(t, lambda, ...) {
  N <- length(t) # length of output
  msd <- rep(0, N)
  # prony coefficients
  rC <- prony.coeff(lambda, ...)
  r <- rC$r
  C <- rC$C
  K <- length(C) # number of modes
  # msd calculations
  if(K > 1) {
    B <- C[2:K]^2/r
    B0 <- sum(B)
    for(ii in 1:(K-1)) {
      msd <- msd + B[ii] * exp(-r[ii]*t)
      #lmsdp <- lmsdp + C[ii+1]^2 / (2*r[ii]) * 2 * r[ii] * exp(-r[ii]*t)
    }
  } else B0 <- 0
  C[1]^2 * t + B0 - msd
  #lmsdp <- lmsdp + C[1]^2
  #lmsdp <- lmsdp / msd * t
  #msd
}
