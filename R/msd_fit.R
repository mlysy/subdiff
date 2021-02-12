#' Sample estimate of mean squared displacement.
#'
#' @template args-Xt
#' @param nlag Number of MSD lags to calculate.
#' @param demean Logical; whether or not to subtract from `Xt` a regression-based estimate of its linear drift.
#' @return Sample MSD vector of length `nlag`.
#'
#' @example examples/Xt_setup.R
#' @example examples/msd_fit.R
#'
#' @export
msd_fit <- function(Xt, nlag, demean = TRUE) {
  ## Xt <- as.matrix(Xt)
  check_Xt(Xt)
  N <- nrow(Xt)
  if(demean) {
    ## mu <- (Xt[N,] - Xt[1,])/(N-1)
    ## Xt <- Xt - ((1:N-1) %o% mu)
    Xt <- matrix(stats::resid(lm(Xt ~ matrix(1:N))), nrow = N)
  }
  if(missing(nlag)) {
    # see [stats::acf()]
    nlag <- floor(10 * (log10(N) - log10(ncol(Xt))))
  }
  nlag <- min(nlag, N-1)
  # output
  msd_empirical(Xt, nlag)
}
