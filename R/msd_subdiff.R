#' Calculate effective subdiffusion time window and parameters from an arbitrary mean square displacement curve.
#'
#' @param msd Vector of sample MSDs.
#' @param tseq Vector of time points at which the MSDs are recorded (see Details).
#' @param rel_tol Relative tolerence for difference between msd and linear fit.
#' @param tmax Logical; if `TRUE` estimate `tmax`, otherwise set it to the last MSD timepoint `tseq[length(tseq)]` (See details).
#' @param log Logical; whether or not the time window is measured on log-scale (See details).
#' @return A vector of length 3 if `tmax = FALSE`, or a vector of length 4 if `tmax = FALSE`.
#' @details
#' Effective subdiffusion time window is defined as the longest time window whose power-law fit (computed by function `msd_ls`) is within a small margin of real msd. It is actually performing the linear regression
#' \preformatted{
#' log(msd[tmin, tmax]) ~ log(D) + alpha * log(tseq[tmin,tmax]).
#' }
#' where `tmin` is the starting point of time window and `tmax` is the end of time window. And the relative error between the true msd and power law approximation is defined as
#' \preformatted{
#' epsilon = max((msd[tmin, tmax] - D * tseq[tmin,tmax]^alpha) / msd[tmin, tmax]).
#' }
#' Since in many experiment `tmax` is out-of-observation, this function allows to set `tmax = tseq[length(tseq)]`, which is the end of experiment. By doing this we only search for `tmin`. In addition, the default value of `log = FALSE` means the length of time window is defined as `tmax - tmin`, but this doesn't give the best MSD fit on the log-log scale, as there are exponentially more points as we move right in the graph, such that the right side of the graph will dominate the fit. Setting `log = TRUE` defines the length of time window as `log(tmax) - log(tmin)` and also applies a log-scaled weight in power fitting (more details in `msd_ls`).
#' This function finds the longest time window whose `epsilon` is smaller than `rel_tol` by applying grid search.
#'
#' @example examples/Xt_setup.R
#' @example examples/msd_fit.R
#' @example examples/msd_subdiff.R
#'
#' @export
msd_subdiff <- function(msd, tseq, rel_tol = 0.05, tmax = FALSE, log = FALSE) {
  ntimes <- length(tseq)
  if(length(msd) != ntimes) stop("msd and tseq have inconsistent dimensions.")
  xx <- tseq # covariate
  ind <- !is.na(msd) # non-NA terms
  if(!tmax) {
    # tmax = N*dt
    Theta <- .tmin_search(msd[ind], xx[ind], rel_tol, log)
    names(Theta) <- c("tmin", "alpha", "D")
  } else {
    #tmax requires computation
    Theta <- .twin_search(msd[ind], xx[ind], rel_tol, log)
    names(Theta) <- c("tmin", "tmax", "alpha", "D")
  }
  Theta
}

.tmin_search <- function(msd, tseq, rel_tol, log) {
  N <- length(tseq)
  for(jj in 2:N-1) {
    msd1 <- msd[jj:N]
    tseq1 <- tseq[jj:N]
    theta <- msd_ls(msd1, tseq1, logw = log)
    sline <- log(theta[2]) + theta[1] * log(tseq1)
    rel_tol1 <- max(abs(exp(sline) - msd1) / msd1)
    if(rel_tol1 < rel_tol) break
  }
  c(tseq[jj], theta)
}

.twin_search <- function(msd, tseq, rel_tol, log) {
  N <- length(tseq)
  stg <- matrix(NA, N, 4)
  for(tmax in N:2) {
    for(jj in 2:tmax-1) {
      msd1 <- msd[jj:tmax]
      tseq1 <- tseq[jj:tmax]
      theta <- msd_ls(msd1, tseq1, logw = log)
      sline <- log(theta[2]) + theta[1] * log(tseq1)
      rel_tol1 <- max(abs(exp(sline) - msd1) / msd1)
      if(rel_tol1 < rel_tol) {
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
  if(all(is.na(dstn))) {
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

