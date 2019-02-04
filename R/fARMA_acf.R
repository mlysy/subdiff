#' ACF of the Autoregressive Moving Average Model with fBM Noise
#'
#' @param alpha TBD
#' @param theta parameter of AR model
#' @param rho parameter of MA model
#' @param dT TBD
#' @param N TBD
#' @return Vector of length N
#' @details
#' \code{Yt = theta Yt-1 + rho Xt + ...}, where \code{Xt} is fBM process with parameter \code{alpha}.
#' @export
farma_acf <- function(alpha, theta, rho, dT, N, m = 30) {
  acf1 <- ma_acf(alpha, c(1-theta-sum(rho), rho), dT, N+m+1)
  acf2 <- ar1_acf(acf1, theta, N, m)
  acf2
}


# ACF of unrestricted fma(q)
ma_acf <- function(alpha, rho, dT, N) {
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
ar1_acf <- function(acf1, theta, N, m) {
  a <- c(1, rep(0, N), theta^(m:1))
  b <- c(acf1, rev(acf1[1:m+1]))
  c <- c(theta^(0:m), rep(0, N+m))
  k1 <- Re(fftw::IFFT(fftw::FFT(b)*fftw::FFT(c)))[1:(N+m+1)]
  acf2 <- Re(fftw::IFFT(fftw::FFT(a)*fftw::FFT(k1)))[1:N]
  acf2
}
