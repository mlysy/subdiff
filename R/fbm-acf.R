#' @title Acf function for the increment of fBM process.
#' @param alpha fBM parameter, range from 0 to 2.
#' @param dT interobservation time.
#' @param N number of samples.
#' @return vector of length N.
#' @export
fbm.acf <- function(alpha, dT, N) {
  if(N == 1) {
    acf <- dT^alpha
  } else {
    acf <- (dT*(0:N))^alpha
    acf <- .5 * (acf[1:N+1] + c(acf[2], acf[1:(N-1)]) - 2*acf[1:N])
  }
  acf
}
