#' @title Acf function for the increment of fBM process.
#' @param alpha fBM parameter, range from 0 to 2.
#' @param dt interobservation time.
#' @param N number of samples.
#' @return vector of length N.
#' @export
fbm.acf <- function(alpha, dt, N){
  tseq <- 1:N * dt
  msd2acf(fbm.msd(tseq = tseq, alpha = alpha))
}