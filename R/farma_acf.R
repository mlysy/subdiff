#' Calculate the ACF for the farma model.
#'
#' Compute the autocovariance of increment process of a fARMA(p,q) model (see **Details**).
#'
#' @param alpha Power law exponent of the fBM process. A scalar between 0 and 2.
#' @param phi A vector of autoregressive (AR) coefficients.
#' @param rho A vector of moving-average (MA) coefficients.
#' @template args-dT
#' @template args-N
#' @param m Number of MA coefficients used for approximating AR coefficients (see **Details**).
#'
#' @template ret-acf
#'
#' @details The fARMA(p,q) model is of following form:
#' \deqn{
#' Y_n = \sum_{i=1}^p \phi_i Y_{n-i} + \sum_{j=0}^q \rho_j X_{n-j}
#' }{
#' Y[n] = \phi_1 Y[n-1] + ... + \phi_p Y[n-p] + \rho_0 X[n] + ... + \rho_q X[n-q]
#' }
#' where residuals \eqn{X[n]} follow the fBM process with parameter `alpha` (See [fbm_acf()]).
#'
#' This function returns the autocovariance of stationary increment process \eqn{\Delta Y[n] = Y[n] - Y[n-1]}.
#'
#' In this function, the autoregressive terms are approximated with `m` moving-average terms, as described in Ling et al (2019).
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
farma_acf <- function(alpha, phi, rho, dT, N, m = 30) {
  if(!phi) {
    acf2 <- .fma_acf(alpha, c(1-sum(rho), rho), dT, N)
  } else {
    acf1 <- .fma_acf(alpha, c(1-sum(phi)-sum(rho), rho), dT, N+m+1)
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
    k1 <- Re(IFFT(FFT(b)*FFT(c)))[1:(N+nlag)]
    acf2 <- Re(IFFT(FFT(a)*FFT(k1)))[1:N]
  }
  acf2
}

# ACF of un-parametrizd far(1)
.ar1_acf <- function(acf1, phi, N, m) {
  a <- c(1, rep(0, N), phi^(m:1))
  b <- c(acf1, rev(acf1[1:m+1]))
  c <- c(phi^(0:m), rep(0, N+m))
  k1 <- Re(IFFT(FFT(b)*FFT(c)))[1:(N+m+1)]
  acf2 <- Re(IFFT(FFT(a)*FFT(k1)))[1:N]
  acf2
}
