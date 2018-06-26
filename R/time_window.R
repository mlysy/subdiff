#' Search for time window in sample mean square displacement curves.
#' 
#' @param msd Vector or matrix of sample MSDs, each column corresponding to a different trajectory.
#' @param tseq Vector of time points at which the MSDs are recorded (see Details).
#' @param error Relative tolerence for difference between msd and linear fit
#' @param tmax Logic, end of experiment time N*dT is used as tmax if \code{FALSE}
#' @param log Logic, length of time window is defined as (tmax - tmin) when \code{FALSE}, and defined as (log(tmax) - log(tmin)) when \code{TRUE}
#' @return Matrix of \code{tmin}, \code{alpha} and \code{D} values, includes column \code{tmax} if tmax is \code{TRUE}
#' @details Time window is defined as the longest successive time period whose linear fitted msd is within a certain tolerence.
#' @export
time_window <- function(msd, tseq, error = 0.05, tmax = FALSE, log = FALSE) {
  yy <- as.matrix(msd)
  npaths <- ncol(yy)
  ntimes <- length(tseq)
  if(nrow(yy) != ntimes) stop("msd and tseq have inconsistent dimensions.")
  xx <- tseq # covariate
  if(!tmax) {
    Theta <- matrix(NA, 3, npaths)
    rownames(Theta) <- c("tmin", "alpha", "D")
    for(ii in 1:npaths) {
      ind <- !is.na(yy[,ii])
      Theta[,ii] <- .tmin_search(yy[ind,ii], xx[ind], error)
    }
  } else {
    Theta <- matrix(NA, 4, npaths)
    rownames(Theta) <- c("tmin", "tmax", "alpha", "D")
    for(ii in 1:npaths) {
      ind <- !is.na(yy[,ii])
      qq <- .twin_search(yy[ind,ii], xx[ind], error, log)
      Theta[,ii] <- qq
    }
    
  }
  if(npaths == 1) Theta <- Theta[,1]
  Theta
}



.tmin_search <- function(msd, tseq, error) {
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

.twin_search <- function(msd, tseq, error, log) {
  N <- length(tseq)
  stg <- matrix(NA, N, 4)
  for(tmax in N:2) {
    for(jj in 2:tmax-1) {
      msd1 <- msd[jj:tmax]
      tseq1 <- tseq[jj:tmax]
      theta <- msd_ls(msd1, tseq1)
      sline <- log(theta[2]) + theta[1] * log(tseq1)
      lmsd1 <- log(msd1)
      error1 <- max(abs(sline - lmsd1) / lmsd1)
      if(error1 < error) {
        stg[tmax,] <- c(tseq[jj], tseq[tmax], theta)
        break
      } else if(jj == tmax-1) {
        stg[tmax,] <- rep(NA, 4)
      }
    }
  }
  if(log) {
    dstn <- log(stg[,2]) - log(stg[,1])
  } else {
    dstn <- stg[,2] - stg[,1]
  }
  # check whether it is full of NA
  if(prod(is.na(dstn))) {
    ans <- rep(NA, 4)
  } else {
    # if there are multiple timewindow of same length
    # we choose the first
    index <-(1:N)[dstn == max(dstn, na.rm = T)]
    index <- index[!is.na(index)]
    if(length(index) > 1) {
      ans <- stg[index[1],]
    } else {
      ans <- stg[index,]
    }
  }
  ans
}

