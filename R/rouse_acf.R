#' ACF of the Rouse-GLE model.
#'
#' Autocorrelation of the Rouse-GLE transient subdiffusion model.
#'
#' @details The Rouse-GLE satisfies the integro-differential equation
#' \deqn{
#' F_t - \int_{-\infty}^t \gamma(t-s) \dot X_s ds = 0,
#' }
#' where \eqn{acf_F(t) = k_B T \gamma(t)} and
#' \deqn{
#' \gamma(t) = 1/K sum_{k=1}^K exp(-\lambda_k * t), \lambda_k = (k/K)^(1/alpha)/tau.
#' }
#'
#' As the temperature \eqn{T} is not supplied, the output is only proportional to the autocorrelation, in the sense that `ACF_X(t) = k_B * T * rouse_acf`.
#'
#' @param alpha Subdiffusion coefficient between 0 and 1.
#' @param tau Shortest timescale of memory kernel.
#' @param K Number of relaxation modes.
#' @param dt Time between obserations.
#' @param N Number of observations.
#' @param ... Additional arguments to pass to [prony_acf()].
#' @return Vector of autocorrelations.
#' @export
rouse_acf <- function(alpha, tau, K, dt, N, ...) {
  lambda <- ((1:K)/K)^(1/alpha)/tau
  ## acf <- prony_acf(lambda, dt, N, ...)
  ## acf * K
  prony_acf(lambda = lambda, nu = 1/K, dt = dt, N = N, ...)
}
