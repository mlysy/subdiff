#' @title Prony MSD
#' @description MSD of a GLE with 0 mass and 0 external potential.
#' @details 
#' $$ vsigma * \int_{-\infty}^t \gamma(t-s) \dot x_s \ud s = k_B T vsigma * F_t
#' where the force noise F_t is a sum of independent OU processes, ie,
#' cov(F_s, F_{t+s}) = gamma(t) = exp(-alpha1 t) + ... + exp(-alphaK t).
#' the resulting process is
#' Y_t = C1 X_1t + ... + CK X_Kt,
#' where X_1t = B_1t and X_it, i > 1 are OU processes,
#' d X_it = -ri X_it dt + d B_it,
#' and all B_it are independent.
#' note that cov(X_i0, X_it) = (2*ri)^{-1} * exp(-ri * t)
#' @param t times at which to calculate the MSD.
#' @param alpha coefficients of the sum of OU processes.
#' @param nSteps number of steps in mode finding (defaults to 100)
#' @param tol tolerance for mode finding (defaults to 0, i.e. always take nSteps)
#' @param r roots that can accelerate computations.
#' @param C coefficients that can accelerate computations, can't pass C without r.
#' @return list containing msd, C, r and lmsdp
#' @export
prony.msd <- function(t, alpha, nSteps = 100, tol = 0, r, C) {
  K <- length(alpha)
  msd <- rep(0, length(t))
  lmsdp <- rep(0, length(t))
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
        
        msd <- msd + C[ii+1]^2 / (2*r[ii]) * 2 * (1 - exp(-r[ii]*t))
        lmsdp <- lmsdp + C[ii+1]^2 / (2*r[ii]) * 2 * r[ii] * exp(-r[ii]*t)
      }
    }
  } else {
    r <- NULL
  }
  C[1] <- sqrt(1/sum(1/alpha))
  msd <- msd + C[1]^2 * t
  lmsdp <- lmsdp + C[1]^2
  lmsdp <- lmsdp / msd * t
  list(msd = msd, r = r, C = C, lmsdp = lmsdp)
}