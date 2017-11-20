#' Normalizing transformations for 1D or 2D variance matrices.
#'
#' @param Sigma Variance matrix on the regular scale.
#' @param rho_scale If \code{TRUE}, the variance matrix is expressed as \code{sigma = sqrt{Sigma[1,1]}} in 1D and \code{(sigma1, sigma2, rho)} in 2D.  Otherwise it is expressed as a \code{1x1} or \code{2x2} matrix.
#' @param D_scale Logical, whether or not the transformation targets the diffusivity constant (see Details).
#' @return Variance matrix on the regular or normalized scale.
#' @details The normalized variance matrix in 1D is \code{lambda = .5 * log(Sigma)}.  In 2D, \code{lambda} is a vector of length 3 with last element \code{lambda[3] = logit(rho/2+1/2)}, where \eqn{\rho} is the correlation.  When \code{D_scale = TRUE}, we have \code{lambda[1:2] = (log(D), log(Sigma[1,1]/Sigma[2,2]))}, where \eqn{D = (\Sigma_{11} + \Sigma_{22})/2}.  Otherwise, the normalizing transformation is \code{lambda[1:2] = (log(Sigma[1,1]), log(Sigma[2,2]))}.
#' @name trans_Sigma
#' @export
trans_Sigma <- function(Sigma, D_scale = TRUE, rho_scale = FALSE) {
  if(length(Sigma) == 1) { # q = 1
    theta <- if(rho_scale) log(Sigma[1]) else .5*log(Sigma[1])
  } else {
    theta <- rep(NA, 3)
    if(rho_scale) {
      theta[1:2] <- Sigma[1:2] # sig1, sig2
      rho <- Sigma[3] # rho
    } else {
      theta[1:2] <- Sigma[c(1,4)] # Sig1, Sig2
      rho <- Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2]) # rho
    }
    if(D_scale) {
      if(rho_scale) theta[1:2] <- theta[1:2]*theta[1:2] # Sig1, Sig2
      theta <- log(c(theta[1:2], .5*(theta[1] + theta[2])))
      theta[1:2] <- c(theta[3], # log(D)
                      theta[1] - theta[2]) # log(Sig1/Sig2)
    } else {
      theta[1:2] <- log(theta[1:2])
      if(!rho_scale) theta[1:2] <- .5*theta[1:2]
    }
    theta[3] <- log((1+rho)/(1-rho)) # logit (rho+1)/2
  }
  theta
}

#' @rdname trans_Sigma
#' @param lambda Variance matrix on the normalized scale (see Details).
#' @export
itrans_Sigma <- function(lambda, D_scale = TRUE, rho_scale = FALSE) {
  if(length(lambda) == 1) {
    theta <- if(rho_scale) exp(lambda) else matrix(exp(2*lambda),1,1) # q = 1
  } else {
    rho <- 2/(1+exp(-lambda[3])) - 1 # 2*ilogit(lambda[3]) - 1
    if(D_scale) {
      Sig <- rep(NA, 2)
      Sig[1] <- exp(lambda[2]) # Sig1/Sig2
      Sig[2] <- 2*exp(lambda[1])/(1 + Sig[1])
      Sig[1] <- Sig[1] * Sig[2]
      if(rho_scale) {
        theta <- c(sqrt(Sig), rho)
      } else {
        theta <- matrix(NA, 2, 2)
        theta[c(1,4)] <- Sig
        theta[1,2] <- theta[2,1] <- rho*sqrt(Sig[1]*Sig[2]) # Sig12
      }
    } else {
      if(rho_scale) {
        theta <- c(exp(lambda[1:2]), rho)
      } else {
        theta <- matrix(NA, 2, 2)
        theta[c(1,4)] <- exp(lambda[1:2]) # sigma
        theta[1,2] <- theta[2,1] <- rho*(theta[1]*theta[4]) # Sig12
        theta[c(1,4)] <- theta[c(1,4)] * theta[c(1,4)] # Sigma
      }
    }
  }
  theta
}

