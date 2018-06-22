#' Search for time window in sample mean square displacement curves.
#' 
#' @param msd Vector or matrix of sample MSDs, each column corresponding to a different trajectory.
#' @param tseq Vector of time points at which the MSDs are recorded (see Details).
#' @param error Relative allowance for difference between msd and linear fit
#' @return A 3-row matrix of \code{tmin}, \code{alpha} and \code{D} values
#' @details Time window is defined as the longest successive time period whose linear fitted msd is within a certain tolerence.
#' Time window is determined by two parameters \code{tmin} and \code{tmax}. 
#' In current version we simply assume that \code{tmax} is always out-of-observation and equals the end of experimental time scale, and focus on the search of \code{tmin}.
#' @export
time_window <- function(msd, tseq, error = 0.05) {
  yy <- as.matrix(msd)
  npaths <- ncol(yy)
  ntimes <- length(tseq)
  if(nrow(yy) != ntimes) stop("msd and tseq have inconsistent dimensions.")
  xx <- tseq # covariate
  Theta <- matrix(NA, 3, npaths)
  rownames(Theta) <- c("tmin", "alpha", "D")
  for(ii in 1:npaths) {
    ind <- !is.na(yy[,ii])
    Theta[,ii] <- .t_search(yy[ind,ii], xx[ind], error)
  }
  if(npaths == 1) Theta <- Theta[,1]
  Theta
}



.t_search <- function(msd, tseq, error) {
  N <- length(tseq)
  for(jj in 2:N-1) {
    msd1 <- msd[jj:N]
    tseq1 <- tseq[jj:N]
    theta <- msd_ls(msd1, tseq1)
    sline <- log(theta[2]) + theta[1] * log(tseq1)
    lmsd1 <- log(msd1)
    error1 <- max(abs(sline - lmsd1) / lmsd1)
    if(error1 < error) break
  }
  c(tseq[jj], theta)
}
