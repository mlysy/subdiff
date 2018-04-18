#' ACF of the Moving Average Model with fBM Noise
#'
#' @param alpha TBD
#' @param rho parameter of MA model
#' @param dT TBD
#' @param N TBD
#' @return Vector of length N
#' @details
#' \code{Yt = (1-rho) Xt + rho Xt-1}, where \code{Xt} is fBM process with parameter \code{alpha}.
#' @export
fma_acf <- function(alpha, rho, dT, N) {
  nlag <- length(rho)
  acf1 <- fbm_acf(alpha, dT, N+nlag)  
  acf2 <- rep(NA, N)
  if(nlag == 1) {
    rho1 <- 1 - rho
    acf2[1] <- (rho1^2 +rho^2) * acf1[1] + 2*rho*rho1*acf1[2]
    acf2[-1] <- (rho1^2 +rho^2) * acf1[-c(1,N+1)] + rho*rho1*(acf1[-c(N,N+1)]+acf1[-c(1,2)])  
  } else {
    for(ii in 1:N) {
      acf2[ii] <- poly_func(rho, ii-1, acf1)
    }
  }
  acf2
}

conv_func <- function(rho, nlag) {
  n <- length(rho)
  sum(rho[1:(n-nlag)] * rho[(1+nlag):n])
}

poly_func <- function(rho, n, acf1) {
  nlag <- length(rho)
  rho1 <- c(1-sum(rho), rho)
  ans <- 0
  for(jj in -nlag:nlag) {
    ans <- ans + conv_func(rho1, abs(jj)) * acf1[abs(n+jj)+1]
  }
  ans
}

