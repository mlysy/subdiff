#' @title Rouse MSD
#' @description mean squared-displacement for the Rouse anomalous diffusion model.
#' @details this corresponds to a GLE with 0 potential and 0 mass with prony series memory kernel 
#' gamma(t) = sum_{n=1}^K 1/K * exp(-alpha_n * t), alpha_n = (n/K)^rho/tau.
#' the output is only proportional to the autocorrelation, in the sense that
#' msd_x(t) = (2 * k_B * T * K) * msd
#' @param t
#' @param rho
#' @param tau
#' @param K
#' @return vector
#' @export
rouse.msd <- function(t, rho, tau, K, ...){
  lambda <- ((1:K)/K)^rho / tau
  rst <- pronyMsd(t, alpha, ...)
  msd <- K * rst$msd
  list(msd = msd, r = rst$r, lmsdp = rst$lmsdp)
}