
#-------------------------------------------------------------------------------

# ok all fitting functions should do the following:
# fit
# coef
# vcov
# resid

# transform parameters to normalized scale:
# gamma = logit(alpha/2)
# nu = log((Sig11 +  Sig22)/2)
# kappa = log(Sig11) - log(Sig22)
# omega = logit((rho+1)/2), rho = Sig12/sqrt(S11*S22)
lsc_trans <- function(alpha, mu, Sigma) {
  q <- length(mu) # dimension check
  theta <- rep(NA, ifelse(q == 1, 3, 6))
  theta[1] <- logit(alpha/2)
  theta[1+1:q] <- mu
  if(q == 1) {
    theta[3] <- log(Sigma[1])/2
  } else if (q == 2) {
    theta[4] <- log((Sigma[1, 1] + Sigma[2, 2])/2)
    theta[5] <- log(Sigma[1, 1]) - log(Sigma[2, 2])
    theta[6] <- logit((Sigma[1, 2]/(sqrt(Sigma[1, 1]*Sigma[2, 2]))+1)/2)
  } else {
    stop("Must have q less than 3.")
  }
  theta
}

lsc_itrans <- function(theta) {
  p <- length(theta)
  q <- ifelse(p == 3, 1, 2)
  alpha <- 2 * ilogit(theta[1])
  Beta <- matrix(theta[1+1:q], 1, q)
  Sigma <- matrix(NA, q, q)
  if(p == 3) { # q = 1
    Sigma[1] <- exp(2 * theta[3])
  } else if (p == 6) { # q = 2
    Sig1 <- exp(theta[5]) # Sig1/Sig2
    Sig2 <- exp(2 * theta[4]) # Sig1+Sig2
    Sigma[2,2] <- Sig2/(1+Sig1) # exp(2*nu)/(1+exp(kappa))
    Sigma[1,1] <- Sig1 * Sigma[2,2]
    rho <- 2*ilogit(theta[6]) - 1
    Sigma[2, 1] <- Sigma[1, 2] <- rho * sqrt(Sigma[1, 1] * Sigma[2, 2])
  } else {
    stop("Must have q less than 3.")
  }
  ans <- list(alpha = alpha, Beta = Beta, Sigma = Sigma)
  ans
}

