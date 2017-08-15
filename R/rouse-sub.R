#' Subdiffusive timescale of Rouse-GLE model.
#'
#' Calculates the range of the subdiffusive timescale along with the effective subdiffusion coefficient and diffusivity constant on this range.
#'
#' @details
#' The Rouse model is a zero-mass, zero-potential GLE with force memory
#' kernel
#' \deqn{
#' acf_F(t) = 1/K sum( exp(- t * (1:K/K)^rho / tau) ).
#' }
#' As \eqn{K} gets large, it exhibits has transient subdiffusion:
#' \deqn{
#' MSD(t) ~ D * t^alpha, t \in tScale
#'          t,           t \notin tScale,
#' }
#' where \eqn{alpha = 1/rho} is the subdiffusion exponent and \eqn{D} is the diffusivity
#' constant.
#' @param alpha Rouse-GLE subdiffusion coefficient.
#' @param K Number of modes in relaxation spectrum.
#' @param tau Shortest timescale of force memory.
#' @param rouse.alpha whether to use the Rouse-GLE's \code{alpha} parameter or estimate \code{alpha} by least-squares.
#' @param ... Additional arguments to pass to \code{prony.coeff}.
#' @return Vector with named elements \code{tmin}, \code{tmax}, \code{alpha}, and \code{D}.
#' @export
rouse.sub <- function(alpha, tau, K, rouse.alpha = FALSE, ...) {
  nPrec <- 5e3 # for fitting D_eff and alpha_eff
  lambda <- ((1:K)/K)^(1/alpha)/tau
  # Rouse coefficients
  rC <- prony.coeff(lambda, ...)
  r <- rC$r
  C <- rC$C
  # subdiffusive timescale
  uScale <- .init.uScale(alpha, tau, K, r, C) # overcover, then decrease
  # inflection method for subdiffusive timescale
  uMin <- optimize(f = function(u) .prony.lmsd(u, r, C)[3],
                   interval = uScale)$min
  uMax <- optimize(f = function(u) -.prony.lmsd(u, r, C)[3],
                   interval = uScale)$min
  # effective subdiffusion and diffusivity coefficients
  if(!rouse.alpha) alpha <- NULL
  aD <- .subdiff.fit(seq(uMin, uMax, len = nPrec), r, C, alpha)
  c(tmin = exp(uMin), tmax = exp(uMax), aD)
}

# initialize the search interval for subdiffusive timescale
# start with smallest/largest timescale of force memory,
# decrease/increase limits by 10x until 1st der of log-log msd > .99
.init.uScale <- function(alpha, tau, K, r, C) {
  uScale <- log(tau * c(1, K)^(1/alpha)) # quick approximation
  #tScale <- tau * c(1, K)^(1/alpha) # quick approximation
  # expand interval until der. of loglog msd ~ 1
  lmsdp <- c(.prony.lmsd(uScale[1], r, C)[2],
             .prony.lmsd(uScale[2], r, C)[2])
  lg10 <- log(10)
  for(ii in 1:100) {
    cutoff <- lmsdp < .99
    if(any(cutoff)) {
      if(cutoff[1]) {
        uScale[1] <- uScale[1] - lg10
        lmsdp[1] <- .prony.lmsd(uScale[1], r, C)[2]
      }
      if(cutoff[2]) {
        uScale[2] <- uScale[2] + lg10
        lmsdp[2] <- .prony.lmsd(uScale[2], r, C)[2]
      }
    } else break
  }
  if(any(lmsdp < .99)) stop("Initial subdiffusive range search failed.")
  uScale
}

# log-log prony-msd and its derivatives.
# u: log-time, a scalar.
.prony.lmsd <- function(u, r, C) {
  K <- length(C)
  t <- exp(u)
  if(K > 1) {
    B <- C[2:K]^2/r
    B0 <- sum(B)
    lmsd <- B * exp(-r * t)
    lmsdp <- r * lmsd
    lmsdp2 <- r * lmsdp
    lmsd <- sum(lmsd)
    lmsdp <- sum(lmsdp)
    lmsdp2 <- sum(lmsdp2)
  } else {
    B0 <- 0
    lmsd <- 0
    lmsdp <- 0
    lmsdp2 <- 0
  }
  A <- C[1]^2
  # calcs on regular scale (i.e., /dt and /dt^2)
  lmsd <- A * t + B0 - lmsd
  lmsdp <- (A + lmsdp)/lmsd
  lmsdp2 <- -(lmsdp2/lmsd + lmsdp^2)
  # calcs on log scale (i.e., /du and /du^2)
  lmsd <- log(lmsd)
  lmsdp <- lmsdp * t
  lmsdp2 <- lmsdp2 * t^2 + lmsdp
  c(x = lmsd, v = lmsdp, a = lmsdp2)
}

# effective subdiffusion and diffusivity constants
# only estimates alpha if it is NULL
.subdiff.fit <- function(useq, r, C, alpha = NULL) {
  msd <- prony.msd(exp(useq), r = r, C = C)*K
  yy <- log(msd)
  xx <- useq
  ybar <- mean(yy)
  xbar <- mean(xx)
  yy <- yy - ybar
  xx <- xx - xbar
  if(is.null(alpha)) {
    alpha <- mean(yy * xx) / mean(xx^2)
  }
  D <- exp(ybar - alpha * xbar)
  c(alpha = alpha, D = D)
}

#--- scratch -------------------------------------------------------------------

# @param tRange how to estimate the subdiffusive timescale.
# @param alphaRouse logical, whether or not to use alpha = 1/rho or the least-squares estimate. Defaults to false.
# @note Either 'force' or 'inflection' method.  When tRange is empty or
# missing, it is the former, which uses tScale = tau * [1, K^(1/alpha)].  When tRange
# is a vector of length 2, the inflection
# method uses the above property of transient subdiffusion, which means
# that on the log-log scale, the slope of the MSD goes from 1 to alpha
# < 1 to 1. Therefore, its acceleration starts at zero, becomes
# negative, then positive, then returns to zero. The two inflection
# points are thus the min and max of d2 log(MSD(t)) / d [log(t)]^2.
# @return list containing tScale(the subdiffusion timescale), alpha(the subdiffusion exponent) and
# D(the diffusivity coefficient)
