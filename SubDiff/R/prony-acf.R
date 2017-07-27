#' @title Prony Acf.
#' @description GLE with 0 mass and 0 external potential.
#' @details
#' process \code{x} is defined as 
#' \code{vsigma * int_{-infty}^t gamma(t-s)  x_s d s = k_B T vsigma * F_t}
#' where the force noise F_t is a sum of independent OU processes, ie,
#' cov(F_s, F_{t+s}) = gamma(t) = exp(-alpha1 t) + ... + exp(-alphaK t).
#' the resulting process is
#' Y_t = C1 X_1t + ... + CK X_Kt,
#' where X_1t = B_1t and X_it, i > 1 are OU processes,
#' d X_it = -ri X_it dt + d B_it,
#' and all B_it are independent.
#' note that cov(X_i0, X_it) = (2*ri)^{-1} * exp(-ri * t).
#' @param lambda coefficients of the sum of OU processes.
#' @param N number of samples.
#' @param dt interobservation time.
#' @param nSteps number of steps in mode finding (defaults to 100).
#' @param tol tolerance for mode finding (defaults to 0, i.e. always take nSteps).
#' @param r roots that can accelerate computations.
#' @param C coefficients that can accelerate computations, can't pass C without r.
#' @note 
#' for backward compatibility, if options is a single number taken to be nSteps and 
#' a warning is issued.
#' @return list with member r: acf(OU_i) propto exp(-ri t) and C: first is weight of the BM, remaining are weights of the OU's.
#' @export
prony.acf <- function(lambda, N, dt, nSteps = 100, tol = 0, r, C) {
  K <- length(lambda) # number of modes
  acf <- rep(0, N)
  hasr <- !missing(r)
  hasC <- !missing(C)
  if(hasC && !hasr) stop("Cannot supply C without r.")
  if(!hasC) C <- rep(0, K)
  # root and coefficient calculations
  if(K > 1) {
    if(!hasr) r <- modePoly(lambda, nSteps, tol)
    if(!hasC) {
      for(ii in 1:(K-1)) {
        tmp <- 1/(lambda-r[ii])
        C[ii+1] <- sqrt(sum(lambda*tmp^2))/r[ii]
        C[ii+1] <- C[ii+1]/(-crossprod(tmp)[1] + sum(tmp)^2)
        
        tmp2 <- exp(-r[ii]*(0:N)*dt)
        tmp3 <- 2*tmp2[1:N] - tmp2[1:N+1] - c(tmp2[2], tmp2[1:(N-1)])
        acf <- acf + C[ii+1]^2/(2*r[ii]) * tmp3
      }
    }
  } else {
    r <- NULL
  }
  if(!hasC) {
    C[1] <- sqrt(1/sum(1/lambda))
    acf[1] <- acf[1] + C[1]^2*dt
  } 
  list(acf = acf, r = r, C = C)
}