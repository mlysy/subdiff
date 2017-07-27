#' @title rouse subdiffusive
#' @description 
#' this function determines the timescale and coefficients of the subdiffusion
#' approximation for the Rouse-GLE model.
#' @details 
#' The Rouse model is a zero-mass, zero-potential GLE with force memory
#' kernel
#' 
#' acf_F(t) = 1/K sum( exp(- t * (1:K/K)^rho / tau) ).
#' 
#' As K gets large, it exhibits has transient subdiffusion:
#' 
#' MSD(t) ~ D * t^alpha, t \in tScale
#'          t,           t \notin tScale,
#'           
#' where alpha = 1/rho is the subdiffusion exponent and D is the diffusivity
#' constant.
#' @param alpha theoretical subdiffusion as K goes to infinity
#' @param K number of modes in relaxation spectrum
#' @param tau relaxation timescale: min and max relaxation times of force F_t are tau * [1, K^(1/alpha)]
#' @param tRange how to estimate the subdiffusive timescale.
#' @param alphaRouse logical, whether or not to use alpha = 1/rho or the least-squares estimate. Defaults to false.
#' @note Either 'force' or 'inflection' method.  When tRange is empty or
#' missing, it is the former, which uses tScale = tau * [1, K^(1/alpha)].  When tRange 
#' is a vector of length 2, the inflection
#' method uses the above property of transient subdiffusion, which means
#' that on the log-log scale, the slope of the MSD goes from 1 to alpha
#' < 1 to 1. Therefore, its acceleration starts at zero, becomes
#' negative, then positive, then returns to zero. The two inflection
#' points are thus the min and max of d2 log(MSD(t)) / d [log(t)]^2.
#' @return list containing tScale(the subdiffusion timescale), alpha(the subdiffusion exponent) and 
#' D(the diffusivity coefficient)
#' @export
rouseSub <- function(rho, tau, K, tRange, alphaRouse = FALSE) {
  nPrec <- 5e3; # should change this to fminsearch...  
  methodForce <- !missing(tRange)
  if(methodForce) {
    # inflection method
    uMin <- log(tRange[1])
    uMax <- log(tRange[2])
    uScale <- c(optimize(f = rouseAcc, interval = c(uMin, uMax))$objective, 
                optimize(f = rouseAcc.sub, interval = c(uMin, uMax))$objective)
    tScale <- exp(uScale)
  } else {
    # force method
    tScale <- tau * c(1, K)^rho
    uScale <- log(tScale)
  }
  
  # get subdiffusion estimates
  tSeq <- exp(seq(from = uScale[1], to = uScale[2], length.out = nPrec))
  msd <- rouseMsd(tseq, rho, tau, K)
  if(alphaRouse) {
    rst <- sDiffEst(tSeq, msd, 1/rho)
  } else {
    rst <- sDiffEst(tSeq, msd)
  }
}

#' nested function rouseAcc
#' 2nd derivative of MSD on log-log scale
rouseAcc <- function(u) {
  t <- exp(u)
  rst <- rouseMsd(t, rho, tau, K)
  acc <- (rst$r + rst$lmsdp * t) * t
  acc
}

#' nestec function rouseAcc.sub
rouseAcc.sub <- function(u) {
  u - rouseAcc(u)
}