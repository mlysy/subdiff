#' Sample estimate of mean squared displacement.
#'
#' @param Xt Vector or matrix of trajectory positions (each column is a different coordinate).  The trajectory is assumed to be sampled at a constant frequency.
#' @param nlag Number of MSD lags to calculate.
#' @param demean Logical; whether or not to remove the mean of \code{Xt} as estimated by linear drift.
#' @return Sample MSD vector of length \code{nlag}.
#' 
#' @example examples/Xt_setup.R
#' @example examples/msd_fit.R
#' 
#' @export
msd_fit <- function(Xt, nlag, demean = TRUE) {
  Xt <- as.matrix(Xt)
  N <- nrow(Xt)
  if(demean) {
    mu <- (Xt[N,] - Xt[1,])/(N-1)
    Xt <- Xt - ((1:N-1) %o% mu)
  }
  if(missing(nlag)) nlag <- floor(10 * (log10(N) - log10(ncol(Xt))))
  nlag <- min(nlag, N-1)
  # output
  .SampleMSD(Xt, nlag)
}
