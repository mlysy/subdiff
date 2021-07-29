#' Calculate the autocorrelation of the fARMA model.
#'
#' Compute the autocorrelation of fARMA(p,q) process increments at equally-spaced time points (see 'Details').
#'
#' @param alpha Subdiffusion exponent of the underlying fBM process. A scalar between 0 and 2.
#' @param phi A vector of `p >= 0` autoregressive (AR) coefficients.
#' @param rho A vector of `q >= 0` moving-average (MA) coefficients (see 'Details').
#' @template args-dt
#' @template args-N
#' @param m Number of MA coefficients used for approximating the ARMA filter (see **Details**).
#'
#' @template ret-acf
#'
#' @details Let `X_n` denote the position of an fBM process at time `t = n * dt`. The fARMA(p,q) process `Y_n` is then defined as
#' ```
#' Y_n = rho_0 X_n + sum_{i=1}^p phi_i Y_{n-i} + sum_{j=1}^q + rho_j X_{n-j}.
#' ```
#' where `rho_0 = 1 - sum(phi) - sum(rho)`.
#'
#' This function returns the autocorrelation of the stationary process `dY_n = Y_{n+1} - Y_n`.  The `ARMA(p,q)` filter is approximated with `m` moving-average terms (see [arma_acf()] for details).
#'
#' @references Ling, Y., Lysy, M., Seim, I., Newby, J.M., Hill, D.B., Cribb, J., and Forest, M.G. "Measurement error correction in particle tracking microrheology" (2019). <https://arxiv.org/abs/1911.06451>.
#'
#' @export
farma_acf <- function(alpha, phi = numeric(), rho = numeric(), dt, N, m = 50) {
  rho0 <- (1 - sum(rho) - sum(phi))
  if(length(rho) != 0) {
    rho <- rho/rho0
  }
  if(length(phi) == 0) {
    # pure MA filter
    facf <- fbm_acf(alpha, dt, N+length(rho))
    acf <- ma_acf(facf, rho = rho)
  } else {
    facf <- fbm_acf(alpha, dt, N+m)
    acf <- arma_acf(facf, phi = phi, rho = rho, m = m)
  }
  # scaling
  rho0^2 * acf
}
