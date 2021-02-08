#' Unconstraining transformation for variance matrices.
#'
#' @param Sigma Variance matrix on the regular scale.
#' @return Variance matrix on the regular or unconstrained scale (see 'Details').
#'
#' @details The unconstraining transformation of a variance matrix is the so-called log-Cholesky decomposition.  Namely, the log-Cholesky decomposition of a variance matrix `Sigma` is a vector `lambda` corresponding to the upper triangular Cholesky factor, of which we take the log of the diagonal and then concatenate the non-zero elements in column-major order.  The exact calculation is given by:
#' ```
#' lambda <- chol(Sigma)
#' diag(lambda) <- log(diag(lambda))
#' lambda <- lambda[upper.tri(lambda,diag=TRUE)]
#' ```
#' The function `trans_Sigma()` converts `Sigma` to `lambda`, whereas `itrans_Sigma()` performs the inverse transformation from `lambda` to `Sigma`.
#' @name trans_Sigma
#' @export
trans_Sigma <- function(Sigma) {
  lambda <- chol(Sigma)
  diag(lambda) <- log(diag(lambda))
  lambda[upper.tri(lambda, diag = TRUE)]
}

#' @rdname trans_Sigma
#' @param lambda Variance matrix on the normalized scale (see 'Details').
#' @export
itrans_Sigma <- function(lambda) {
  # determine size of matrix
  ndim <- (-1 + sqrt(1 + 8*length(lambda)))/2
  # construct upper Cholesky factor
  U <- matrix(0, ndim,ndim)
  U[upper.tri(U, diag = TRUE)] <- lambda
  diag(U) <- exp(diag(U))
  # construct variance matrix
  crossprod(U)
}

#--- scratch -------------------------------------------------------------------

## #' Unconstraining transformations for 1D or 2D variance matrices.
## #'
## #' @param Sigma Variance matrix on the regular scale.
## #' @param rho_scale If `TRUE`, the variance matrix is expressed as \code{sigma = sqrt{Sigma[1,1]}} in 1D and `(sigma1, sigma2, rho)` in 2D.  Otherwise it is expressed as a `1x1` or `2x2` matrix.
## #' @param D_scale Logical, whether or not the transformation targets the diffusivity constant (see Details).
## #' @return Variance matrix on the regular or normalized scale.
## #' @details The unconstrained variance matrix is defined as follows:
## #' \itemize{
## #'   \item In 1D, it is `lambda = .5 * log(Sigma)`.
## #'   \item In 2D, `lambda` is a vector of length 3 with last element `lambda[3] = logit(rho/2+1/2)`, where \eqn{\rho = \Sigma_{12}/} is the correlation.  When `D_scale = TRUE`, we have `lambda[1:2] = (log(D), log(Sigma[1,1]/Sigma[2,2]))`, where \eqn{D = (\Sigma_{11} + \Sigma_{22})/2}.  Otherwise, the unconstraining transformation is `lambda[1:2] = (log(Sigma[1,1]), log(Sigma[2,2]))`.
## #' }
## #' @name trans_Sigma
## #' @export
## trans_Sigma <- function(Sigma, D_scale = TRUE, rho_scale = FALSE) {
##   if(length(Sigma) == 1) { # q = 1
##     theta <- if(rho_scale) log(Sigma[1]) else .5*log(Sigma[1])
##   } else {
##     theta <- rep(NA, 3)
##     if(rho_scale) {
##       theta[1:2] <- Sigma[1:2] # sig1, sig2
##       rho <- Sigma[3] # rho
##     } else {
##       theta[1:2] <- Sigma[c(1,4)] # Sig1, Sig2
##       rho <- Sigma[1,2]/sqrt(Sigma[1,1]*Sigma[2,2]) # rho
##     }
##     if(D_scale) {
##       if(rho_scale) theta[1:2] <- theta[1:2]*theta[1:2] # Sig1, Sig2
##       theta <- log(c(theta[1:2], .5*(theta[1] + theta[2])))
##       theta[1:2] <- c(theta[3], # log(D)
##                       theta[1] - theta[2]) # log(Sig1/Sig2)
##     } else {
##       theta[1:2] <- log(theta[1:2])
##       if(!rho_scale) theta[1:2] <- .5*theta[1:2]
##     }
##     theta[3] <- log((1+rho)/(1-rho)) # logit (rho+1)/2
##   }
##   theta
## }

## #' @rdname trans_Sigma
## #' @param lambda Variance matrix on the normalized scale (see 'Details').
## #' @export
## itrans_Sigma <- function(lambda, D_scale = TRUE, rho_scale = FALSE) {
##   if(length(lambda) == 1) {
##     theta <- if(rho_scale) exp(lambda) else matrix(exp(2*lambda),1,1) # q = 1
##   } else {
##     rho <- 2/(1+exp(-lambda[3])) - 1 # 2*ilogit(lambda[3]) - 1
##     if(D_scale) {
##       Sig <- rep(NA, 2)
##       Sig[1] <- exp(lambda[2]) # Sig1/Sig2
##       Sig[2] <- 2*exp(lambda[1])/(1 + Sig[1])
##       Sig[1] <- Sig[1] * Sig[2]
##       if(rho_scale) {
##         theta <- c(sqrt(Sig), rho)
##       } else {
##         theta <- matrix(NA, 2, 2)
##         theta[c(1,4)] <- Sig
##         theta[1,2] <- theta[2,1] <- rho*sqrt(Sig[1]*Sig[2]) # Sig12
##       }
##     } else {
##       if(rho_scale) {
##         theta <- c(exp(lambda[1:2]), rho)
##       } else {
##         theta <- matrix(NA, 2, 2)
##         theta[c(1,4)] <- exp(lambda[1:2]) # sigma
##         theta[1,2] <- theta[2,1] <- rho*(theta[1]*theta[4]) # Sig12
##         theta[c(1,4)] <- theta[c(1,4)] * theta[c(1,4)] # Sigma
##       }
##     }
##   }
##   theta
## }

