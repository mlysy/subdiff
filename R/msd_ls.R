#' Least-squares power-law fit to sample mean square displacement curves.
#'
#' @param msd Vector or matrix of sample MSDs, each column corresponding to a different trajectory.
#' @param tseq Vector of time points at which the MSDs are recorded (see Details).
#' @param pooled Logical; whether to calculate separate regressions for each MSD or pool all of them together.
#' @param logw Logical; whether or not regression should be log-weighted (see Details).
#' @return A vector of length 2 if \code{pooled = TRUE}, or a 2-row matrix of \code{alpha} and \code{D} values if \code{pooled = FALSE}.
#' @details A power-law fit to the empirical MSD is calculated by performing the linear regression
#' \preformatted{
#' log(msd) ~ log(D) + alpha * log(tseq).
#' }
#' From a curve fitting perspective, we seek to estimate the "best" power-law describing the MSD, i.e., if MSD(t) is the true MSD at time t, then the best power law minimizes
#' \preformatted{
#' integral_0^Inf | log MSD(t) - log(D) - alpha * log(t) |^2 w(t) dt,
#' }
#' where \code{w(t) > 0} is a weight function.  The default value of \code{logw = FALSE} uses \code{w(t) = 1}, but this doesn't give the best MSD fit on the log-log scale, as there are exponentially more points as we move right in the graph, such that the right side of the graph will dominate the fit.  Setting \code{logw = TRUE} uses uniform weighting on the log-scale, which corresponds to \code{w(t) = 1/t}.
#'
#' For computational efficiency, the MSDs are expected to all be sampled at the same time points.  If this is not the case, missing time points should beindicated by \code{NA} in the \code{msd} matrix.
#' 
#' @example examples/Xt_setup.R
#' @example examples/msd_fit.R
#' @example examples/msd_ls.R
#' 
#' @export
msd_ls <- function(msd, tseq, pooled = TRUE, logw = TRUE) {
  yy <- log(as.matrix(msd))
  npaths <- ncol(yy)
  ntimes <- length(tseq)
  if(nrow(yy) != ntimes) stop("msd and tseq have inconsistent dimensions.")
  xx <- log(tseq) # covariate
  # regression weights
  ww <- if(logw) 1/tseq else rep(1, length(tseq))
  if(pooled && (npaths > 1)) {
    # combined weights/response for pooled estimate
    ww <- ww * apply(!is.na(yy), 1, sum)
    yy <- as.matrix(rowMeans(yy, na.rm = TRUE))
    npaths <- 1
  }
  Theta <- matrix(NA, 2, npaths)
  rownames(Theta) <- c("alpha", "D")
  for(ii in 1:npaths) {
    ind <- !is.na(yy[,ii])
    Theta[,ii] <- .weighted_ls(yy[ind,ii], xx[ind], ww[ind])
  }
  if(npaths == 1) Theta <- Theta[,1]
  Theta
}

# efficient weighted regression in alpha/D estimation context.
# will eventually be converted to C++
.weighted_ls <- function(yy, xx, ww) {
  ww <- ww/sum(ww) # normalize weights
  ybar <- sum(ww * yy)
  xbar <- sum(ww * xx)
  yy <- yy - ybar
  xx <- xx - xbar
  alpha <- sum(ww * yy * xx) / sum(ww * xx^2)
  c(alpha, exp(ybar - alpha * xbar))
}
