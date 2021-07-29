#' Coefficients of position-process representation of the Prony-GLE.
#'
#' @template args-lambda
#' @template args-nu
#' @param r Optional vector of pre-computed OU mean-reversion parameters.  Avoids calling the mode-finding routine.
#' @param C Optional vector of pre-computed OU scale coefficients.  Can't pass `C` without `r`.  When both `r` and `C` are supplied the function does nothign.  This is mainly for the convenience of calling `prony_coeff()` from other functions.
#' @param nsteps Number of steps in mode-finding golden search algorithm.
#' @param tol Relative tolerance in mode-finding golden search algorithm.  For numerical stability, default is 0 such that `nsteps` are always used.
#' @return A list with elements `r` and `C`, containing the `length(lambda)-1` and `length(lambda)` vectors of mean-reversion parameters and scale factors for the BM + OU representation of the Prony-GLE.
#' @template details-prony_gle
#' @note Since temperature is not provide, `C` is only proportional to the true vector `C_true`, such that `C_true = sqrt(k_B T) * C`.
#' @export
prony_coeff <- function(lambda, nu, r, C, nsteps = 100, tol = 0) {
  hasr <- !missing(r)
  hasC <- !missing(C)
  if(hasC) {
    if(!hasr) {
      stop("Cannot supply C without r.")
    } else {
      # do nothing
      return(list(r = r, C = C))
    }
  }
  K <- length(lambda) # number of modes
  C <- rep(0, K)
  # OU components
  if(K > 1) {
    if(!hasr) r <- mode_poly(lambda, nsteps, tol)
    for(ii in 1:(K-1)) {
      tmp <- 1/(lambda-r[ii])
      C[ii+1] <- sqrt(sum(lambda*tmp^2))/r[ii]
      C[ii+1] <- C[ii+1]/(-crossprod(tmp)[1] + sum(tmp)^2)
    }
  } else r <- NULL
  C[1] <- sqrt(1/sum(1/lambda)) # BM component
  list(r = r, C = sqrt(2/nu) * C)
}
