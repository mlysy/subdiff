#' ACF of the Rouse-GLE model.
#'
#' Autocorrelation of the Rouse-GLE transient subdiffusion model.
#' @details The Rouse-GLE satisfies the integro-differential equation
#' \deqn{
#' F_t - \int_{-\infty}^t \gamma(t-s) \dot X_s ds = 0,
#' }
#' where \eqn{acf_F(t) = k_B T \gamma(t)} and
#' \deqn{
#' \gamma(t) = 1/K sum_{n=1}^K exp(-\lambda_n * t), \lambda_n = (n/K)^(1/alpha)/tau.
#' }
#'
#' As the temperature \eqn{T} is not supplied, the output is only proportional to the autocorrelation, in the sense that \eqn{acf_x(t) = (2 * k_B * T * K) * acf}.
#' @param alpha Subdiffusion coefficient between 0 and 1.
#' @param tau Shortest timescale of memory kernel.
#' @param K Number of relaxation modes.
#' @param N Number of observations.
#' @param dt Time between obserations.
#' @param ... Additional arguments to pass to \code{\link{prony_acf}}.
#' @return Vector of autocorrelations.
#' @export
rouse_acf <- function(alpha, tau, K, N, dt, ...) {
  lambda <- ((1:K)/K)^(1/alpha)/tau
  acf <- prony_acf(lambda, N, dt, ...)
  acf * K
}
