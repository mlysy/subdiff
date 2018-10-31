#' ACF of the Autoregressive Moving Average Model with fBM Noise
#'
#' @param alpha TBD
#' @param gamma parameter of AR model
#' @param rho parameter of MA model
#' @param dT TBD
#' @param N TBD
#' @return Vector of length N
#' @details
#' \code{Yt = gamma Yt-1 + rho Xt + (1-gamma-rho) Xt-1}, where \code{Xt} is fBM process with parameter \code{alpha}.
#' @export
farma_acf <- function(alpha, gamma, rho, dT, N, p = 50) {
  p1 <- (1-rho-gamma)/(1-gamma)
  Z_acf <- fma_acf(alpha, p1, dT, N+p+1) * (1-gamma)^2
  a <- c(1, rep(0, N), gamma^(p:1))
  b <- c(Z_acf, rev(Z_acf[1:p+1]))
  c <- c(gamma^(0:p), rep(0, N+p))
  k1 <- Re(fftw::IFFT(fftw::FFT(b)*fftw::FFT(c)))[1:(N+p+1)]
  acf2 <- Re(fftw::IFFT(fftw::FFT(a)*fftw::FFT(k1)))[1:N]
  acf2
}

# Acf of fARMA model with parameters theta and rho
# Y_n = theta Y_n-1 + rho X_n + (1-rho-theta) X_n-1, |theta| < 1
# Y_n = theta Y_n-1 + Zn, 
# Z_n = (1-theta) * [rho/(1-theta) X_n + (1-rho-theta)/(1-theta) X_n-1] has acf
# fma_acf(alpha, (1-rho-theta)/(1-theta), dT, N) * (1-theta)^2

# circulant <-function(x) {
#   n <- length(x)
#   ans <- suppressWarnings(
#     matrix(x[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n))
#   t(ans)
# }
# 
# N <- 6
# p <- 4
# 
# theta <- .8
# a <- c(1, rep(0, N), theta^(p:1))
# A1 <- circulant(a)
# A <- A1[1:(N+1),]
# 
# Z_acf <- as.numeric(0:(N+p))
# b <- c(Z_acf, rev(Z_acf[1:p+1]))
# B1 <- circulant(b)
# B <- B1[1:(N+p+1), 0:p+1]
# 
# c <- theta^(0:p)
# c0 <- c(c, rep(0, N+p))
# 
# ans1 <- c(A %*% B %*% as.matrix(c))
# 
# k1 <- Re(IFFT(FFT(b)*FFT(c0)))[1:(N+p+1)]
# ans2 <- Re(IFFT(FFT(a)*FFT(k1)))[1:(N+1)]
# 
# range(ans1 - ans2)