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
  acf1 <- fbm_acf(alpha, dT, N+1)
  acf2 <- rep(NA, N)
  rho1 <- 1 - rho
  acf2[1] <- (rho1^2 +rho^2) * acf1[1] + 2*rho*rho1*acf1[2]
  acf2[-1] <- (rho1^2 +rho^2) * acf1[-c(1,N+1)] + rho*rho1*(acf1[-c(N,N+1)]+acf1[-c(1,2)])
  acf2
}


# # for ma, Yt = (1-rho) Xt + rho Xt-1
# mat_func <- function(rho, N) {
#   mat <- matrix(0, N, N+1)
#   for(ii in 1:N) {
#     mat[ii, ii] <- 1 - rho
#     mat[ii, ii+1] <- rho
#   }
#   mat
# }
#
# # test
# rho <- .2
# N <- 5
# m1 <- mat_func(rho, N)
# acf1 <- rnorm(N+1)
# m1 %*% toeplitz(acf1) %*% t(m1)
#
# #
#
# range(toeplitz(ar_acf(acf1, rho)) - m1 %*% toeplitz(acf1) %*% t(m1))

