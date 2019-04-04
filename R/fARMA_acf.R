#' Calculate the autocorrelation function for the farma(p,q) model.
#'
#' Compute the autocorrelation of increments of farma(p,q) model (see \strong{Details}).
#' 
#' @param alpha Power law exponent of fBM process. A scalar between 0 and 2.
#' @param phi A vector of AR coefficients.
#' @param rho A vector of MA coefficients.
#' @template args-dT
#' @template args-N
#' @param m Number of MA coefficients used for approximating AR coefficients (see \strong{Details}).
#' 
#' @template ret-acf
#' 
#' @details The farma(p,q) model is of following form:
#' \deqn{
#' Y_t = \sum_{i=1}^p \phi_i Y_{t-i} + \sum_{j=0}^q \rho_j X_{t-j}, \rho_0 = 1-\sum_{i=1}^p \phi_i-\sum_{j=1}^q \rho_j,
#' }{
#' 
#' }
#' where residuals \eqn{X_t} follows the fBM model with parameter \code{alpha}.
#' Auto-regressive terms are approximated with \code{m} terms of Moving-average terms.
#'
#' @examples 
#' acf1 <- farma_acf(alpha = 0.8, phi = 0.1, rho = c(0.1, 0.05), dT = 1/60, N = 1800)
#' 
#' @export
farma_acf <- function(alpha, phi, rho, dT, N, m = 30) {
  if(!phi) {
    acf2 <- .fma_acf(alpha, c(1-sum(rho), rho), dT, N)
  } else {
    acf1 <- .fma_acf(alpha, c(1-phi-sum(rho), rho), dT, N+m+1)
    acf2 <- .ar1_acf(acf1, phi, N, m)
  }
  acf2
}

# ACF of unparametrized moving-average model with fBM noises.
.fma_acf <- function(alpha, rho, dT, N) {
  nlag <- length(rho)
  acf1 <- fbm_acf(alpha, dT, N+nlag)  
  if(nlag == 2) {
    acf2 <- sum(rho^2)*acf1[1:N] + (rho[1]*rho[2])*(acf1[c(2,2:N-1)]+acf1[1:N+1])
  } else if(nlag == 3) {
    acf2 <- sum(rho^2)*acf1[1:N] + 
      (rho[1]*rho[2]+rho[2]*rho[3])*(acf1[c(2,2:N-1)]+acf1[1:N+1]) + 
      (rho[1]*rho[3])*(acf1[c(3,2,3:N-2)]+acf1[1:N+2])
  } else {
    a <- c(rho[1], rep(0, N), rev(rho[-1]))
    b <- c(acf1, rev(acf1[2:nlag]))
    c <- c(rho, rep(0, N+nlag-1))
    k1 <- Re(fftw::IFFT(fftw::FFT(b)*fftw::FFT(c)))[1:(N+nlag)]
    acf2 <- Re(fftw::IFFT(fftw::FFT(a)*fftw::FFT(k1)))[1:N]
  }
  acf2
}

# ACF of un-parametrizd far(1)
.ar1_acf <- function(acf1, phi, N, m) {
  a <- c(1, rep(0, N), phi^(m:1))
  b <- c(acf1, rev(acf1[1:m+1]))
  c <- c(phi^(0:m), rep(0, N+m))
  k1 <- Re(fftw::IFFT(fftw::FFT(b)*fftw::FFT(c)))[1:(N+m+1)]
  acf2 <- Re(fftw::IFFT(fftw::FFT(a)*fftw::FFT(k1)))[1:N]
  acf2
}
