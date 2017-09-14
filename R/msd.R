#' @title Mean Square Displacement
#' @param Xt trajectory
#' @param lag order of MSD of interested, default 20
#' @return vector of length \code{lag}
#' @export
msd <- function(Xt, lag = 20) {
  N <- length(Xt)
  ans <- rep(NA, lag)
  for(ii in 1:lag) {
    ans[ii] <- mean(diff(Xt, lag = ii)^2)
  }
  ans
}