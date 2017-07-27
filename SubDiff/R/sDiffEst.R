#' @title sDiffEst
#' @description this function calculates subdiffusion power law 
#' MSD(t) = D t^alpha
#'  from a discretized MSD least-squares estimate.
#' @details 
#'  [alpha, D] = sDiffEst(tSeq, msd, wgt) calculates these coefficients
#'  from the regression
#'  
#'  'log(msd) ~  log(D) + alpha tSeq + wgt^(-1/2) eps.'
#'  
#'  When tSeq is of length 2 it denotes a range.
#'  alpha if provided is fixed.
#'  wgt if not provided is ones.
#' @param tSeq
#' @param msd
#' @param alpha
#' @wgt
#' @return 
#' @export
sDiffEst <- function(tSeq, msd, alpha, wgt) {
  if(length(tSeq) == 2) {
    tScal <- tSeq
    tSeq <- seq(from = log(tSeq[1]), to = log(tSeq[2]), length.out = length(msd))
  } else {
    tScale <- tSeq
    tSeq <- log(tSeq)
  }
  noW <- missing(wgt)
  if(!noW) {
    wgt <- sqrt(wgt)
  }
  calcAlpha <- missing(alpha)
  
  # regression vectors
  yy <- log(msd)
  xx <- tSeq
  if(!noW) {
    yy <- yy * wgt
    xx <- xx * wgt
  }
  ybar <- mean(yy)
  xbar <- mean(xx)
  yy <- yy - ybar
  xx <- xx - xbar
  
  # coefficients
  if(calcAlpha) {
    alpha = sum(yy * xx) / sum(xx^2)
  }
  D <- exp(ybar - alpha * xbar)
}
