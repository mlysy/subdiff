#' ACF of the Autoregressive Model
#' @param alpha
#' @param rho parameter of MA model
#' @param dT
#' @param N
#' @return vector of length N
#' @details 
#' Yt = (1-rho) Yt-1 + rho Xt, where Xt is fBM process with parameter alpha
#' @export
ar_acf <- function(alpha, rho, dT, N, nlim = 1e2) {
  if(rho > 0) {
    rho1 <- 1-rho
  } else {
    rho1 <- 1/(1-rho)
  }
  acf1 <- fbm_acf(alpha, dT, N+nlim)
  acf2 <- rep(NA, N)
  for(ii in 1:N) {
    acf2[ii] <- .ar_terms(acf1, rho1, ii-1, nlim)
  }
  acf2 * rho^2
}

#' nested function
.ar_terms <- function(acf1, rho1, k, nlim) {
  mat.temp <- matrix(NA, nlim, nlim)
  for(ii in 1:nlim-1) {
    for(jj in 1:nlim-1) {
      mat.temp[ii+1, jj+1] <- rho1^(ii+jj) * acf1[abs(k+ii-jj)+1]
    }
  }
  sum(mat.temp)
}
