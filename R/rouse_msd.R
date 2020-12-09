#' MSD of the Rouse-GLE model.
#'
#' Mean Squared Displacement of the Rouse-GLE transient subdiffusion model.
#'
#' @details The Rouse-GLE satisfies the integro-differential equation
#' \deqn{
#' F_t - \int_{-\infty}^t \gamma(t-s) \dot X_s ds = 0,
#' }
#' where \eqn{acf_F(t) = k_B T \gamma(t)} and
#' \deqn{
#' \gamma(t) = 1/K sum_{n=1}^K exp(-\lambda_n * t), \lambda_n = (n/K)^(1/alpha)/tau.
#' }
#'
#' As the temperature \eqn{T} is not supplied, the output is only proportional to the MSD, in the sense that \eqn{MSD_x(t) = (2 * k_B * T * K) * msd}.
#' @param t Vector of timepoints at which to calculate the MSD.
#' @param alpha Subdiffusion coefficient between 0 and 1.
#' @param tau Shortest timescale of memory kernel.
#' @param K Number of relaxation modes.
#' @param ... Additional arguments to pass to [prony_msd()].
#' @return Vector of mean square displacements.
#' @export
rouse_msd <- function(t, alpha, tau, K, ...) {
  lambda <- ((1:K)/K)^(1/alpha) / tau
  ## msd <- prony_msd(t, lambda, ...)
  ## msd * K
  prony_msd(t = t, lambda = lambda, nu = 1/K, ...)
}
