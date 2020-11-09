#' Calculate the ACF for the farma model.
#'
#' Compute the autocovariance of increment process of a fARMA(p,q) model (see **Details**).
#'
#' @param alpha Power law exponent of the fBM process. A scalar between 0 and 2.
#' @param phi A vector of `p >= 0` autoregressive (AR) coefficients.
#' @param rho A vector of `q >= 0` moving-average (MA) coefficients (see 'Details').
#' @template args-dT
#' @template args-N
#' @param m Number of MA coefficients used for approximating the ARMA filter (see **Details**).
#'
#' @template ret-acf
#'
#' @details The fARMA(p,q) model is of following form:
#' \deqn{
#' Y_n = \sum_{i=1}^p \phi_i Y_{n-i} + \sum_{j=0}^q \rho_j X_{n-j}
#' }{
#' Y[n] = \phi_1 Y[n-1] + ... + \phi_p Y[n-p] + \rho_0 X[n] + ... + \rho_q X[n-q]
#' }
#' where residuals \eqn{X[n]} follow the fBM process with parameter `alpha` (see [fbm_acf()]), and `rho_0 = 1 - sum(phi) - sum(rho)`.
#'
#' This function returns the autocovariance of stationary increment process \eqn{\Delta Y[n] = Y[n] - Y[n-1]}.  The `ARMA(p,q)` filter is approximated with `m` moving-average terms using the method of Brockwell & Davis (1991).
#'
#' @references Ling, Y., Lysy, M., Seim, I., Newby, J.M., Hill, D.B., Cribb, J., and Forest, M.G. "Measurement error correction in particle tracking microrheology." (2019). <https://arxiv.org/abs/1911.06451>.
#'
#' @example examples/acf_setup.R
#' @example examples/farma_acf.R
#'
#' @importFrom fftw FFT
#' @importFrom fftw IFFT
#'
#' @export
farma_acf <- function(alpha, phi = numeric(), rho = numeric(), dT, N, m = 30) {
  rho0 <- (1 - sum(rho) - sum(phi))
  if(length(rho) != 0) {
    rho <- rho/rho0
  }
  if(length(phi) == 0) {
    # pure MA filter
    facf <- fbm_acf(alpha, dT, N+length(rho))
    acf <- ma_acf(facf, rho = rho)
  } else {
    facf <- fbm_acf(alpha, dT, N+m)
    acf <- arma_acf(facf, phi = phi, rho = rho, m = m)
  }
  # scaling
  rho0^2 * acf
}
