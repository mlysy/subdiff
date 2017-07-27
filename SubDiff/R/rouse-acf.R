#' @title Rouse ACF
#' @description 
#' autocorrelation for the Rouse anomalous diffusion model. This corresponds to a GLE with 0 potential and 0 mass with prony series memory kernel.
#' @details 
#' gamma(t) = sum_{n=1}^K 1/K * exp(-alpha_n * t), alpha_n = (n/K)^rho/tau.
#' the output is only proportional to the autocorrelation, in the sense that
#' acf_x(t) = (2 * k_B * T * K) * acf.
#' @param alpha why use 1/alpha in this function?
#' @param tau
#' @param K
#' @param N
#' @param dt
#' @export
rouse.acf <- function(alpha, tau, K, N, dt, ...) {
  lambda <- ((1:K)/K)^(1/alpha)/tau
  acf <- prony.acf(lambda, N, dt, ...)$acf
  acf * K
}