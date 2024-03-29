#' Calculate the autocorrelation of an ARMA(p,q) filter applied to a stationary process.
#'
#' @param acf Vector of length `N + m` of autocorrelations of the original stationary process.
#' @param phi Vector of `p >= 0` coefficients defining the autoregressive part of the filter (see 'Details').
#' @param rho Vector of `q >= 0` coefficients defining the moving-average part of the filter (see 'Details').
#' @param m Order of the MA approximation to the ARMA(p,q) filter.
#' @return A vector of length `N` giving the autocorrelation of the filtered process.
#' @details Let `X_n` denote the observation of the original process at time `n`, and `Y_n` denote the corresponding observation of the filtered process.  The AR(p,q) filter model is defined as
#' ```
#' Y_n = X_n + sum_{i=1}^p phi_i Y_{n-i} + sum_{j=1}^q rho_j X_{n-j}.
#' ```
#' The autocorrelation of `Y_n` can be computed exactly from that of `X_n` for pure moving-average filters with `p = 0`.  The general ARMA(p,q) filter is first approximated by a moving-average process of order `m` of the form
#' ```
#' Y[n] = X[n] + sum_{j=1}^m psi[j] X[n-j],
#' ```
#' where the coefficients `psi` are determined using the method of Brockwell & Davis (1991) implemented in [stats::ARMAtoMA()].
#'
#' @importFrom fftw FFT IFFT
#' @importFrom stats ARMAtoMA
#'
#' @references Brockwell, P.J. and Davis, R.A. "Time Series: Theory and Methods" (1991).  Springer, New York.  <https://doi.org/10.1007/978-1-4899-0004-3>.
#'
#' @export
arma_acf <- function(acf, phi = numeric(), rho = numeric(), m) {
  ## if(length(rho) == 0) stop("rho must have non-zero length.")
  ## if(length(phi) == 0) {
  ##   psi <- rho
  ## } else {
  ##   psi <- arma2ma(phi = phi, rho = rho, q = m)
  ## }
  psi <- ARMAtoMA(ar = phi, ma = rho, lag.max = m)
  if(length(acf) <= m) stop("acf must be of length greater than m.")
  ma_acf(acf, rho = psi)
}

#' Calculate the autocorrelation of an `MA(q)` filter applied to a stationary process.
#'
#' @param acf Vector of length `N + q` of autocorrelations of the original stationary process.
#' @param rho Vector of `q >= 0` coefficients of the `MA(q)` filter (see 'Details').
#' @return A vector of length `N` giving the autocorrelation of the filtered process.
#'
#' @details Let `X_n` denote the observation of the original process at time `n`, and `Y_n` denote the corresponding observation of the filtered process.  The MA(q) filter model is defined as
#' ```
#' Y_n = X_n + sum_{j=1}^q rho_j X_{n-j}.
#' ```
#' @export
ma_acf <- function(acf, rho = numeric()) {
  q <- length(rho)
  N <- length(acf) - q
  if(q == 0) {
    facf <- acf[1:N]
  } else if(q == 1) {
    facf <- (1 + rho^2)*acf[1:N] + rho[1]*(acf[c(2,2:N-1)]+acf[1:N+1])
  } else if(q == 2) {
    facf <- (1 + sum(rho^2))*acf[1:N] +
      (rho[1]+rho[1]*rho[2])*(acf[c(2,2:N-1)]+acf[1:N+1]) +
      rho[2]*(acf[c(3,2,3:N-2)]+acf[1:N+2])
  } else {
    aa <- c(1, rep(0, N), rev(rho))
    bb <- c(acf, rev(acf[1+1:q]))
    cc <- c(1, rho, rep(0, N+q-1))
    dd <- Re(IFFT(FFT(bb)*FFT(cc)))[1:(N+q+1)]
    facf <- Re(IFFT(FFT(aa)*FFT(dd)))[1:N]
  }
  facf
}
