#' @title transform parameters into standard form
#' @description transfrorm {alpha, Beta, Sigma} into {logit(alpha/2), Beta, SLR(Sigma)}, where 
#' SLR(Sigma) = {D = (var1 + var2)/2, D2 = (var1 - var2), logit((rho+1)/2)}
#' @param alpha fBM parameter
#' @param Beta linear drift parameter
#' @param Sigma column-wise variance
#' @param trans logic, transform Matrix Normal parameters into standard form, otherwise just display
#' @return vector
#' @export
transFunc <- function(alpha, Beta, Sigma, trans) {
  # dimension check
  q <- length(Beta)
  if(q == 1) {
    theta <- rep(NA, 3)
    theta[1] <- alpha
    theta[2] <- Beta[1]
    theta[3] <- Sigma[1]
    if(trans) {
      theta[1] <- logit(alpha/2)
      theta[3] <- log(Sigma[1])/2
    }
  } else if (q == 2) {
    theta <- rep(NA, 6)
    theta[1] <- alpha
    theta[2:3] <- Beta[1:2]
    theta[4] <- (Sigma[1, 1] + Sigma[2, 2])/2
    theta[5] <- (Sigma[1, 1] - Sigma[2, 2])/2
    theta[6] <- Sigma[1, 2]
    if(trans) {
      theta[1] <- logit(theta[1]/2)
      theta[4] <- log(Sigma[1, 1])/2
      theta[5] <- log(Sigma[2, 2])/2
      theta[6] <- logit((Sigma[1, 2]/(sqrt(Sigma[1, 1]*Sigma[2, 2]))+1)/2)
    }
  } else {
    stop("So far we can only deal with dimension < 3")
  }
  theta
}


#' @title Inverse function of transFunc
#' @param theta transformed parameters
#' @param trans logic, true means \code{theta} is transformed
#' @return list containing {alpha, Beta, Sigma}
#' @export
itransFunc <- function(theta, trans) {
  p <- length(theta)
  if(p == 3) {
    # q <- 1
    alpha <- theta[1]
    Beta <- matrix(theta[2], 1, 1)
    Sigma <- matrix(theta[3], 1, 1)
    if(trans) {
      alpha <- 2 * ilogit(alpha)
      Sigma[1] <- exp(2 * theta[3])
    } 
  } else if (p == 6) {
    # q <- 2
    alpha <- theta[1]
    Beta <- matrix(theta[2:3], 1, 2)
    Sig1 <- theta[4] + theta[5]
    Sig2 <- theta[4] - theta[5]
    Sigma <- matrix(c(Sig1, theta[6], theta[6], Sig2), 2, 2)
    if(trans) {
      alpha <- 2 * ilogit(alpha)
      Sigma[1, 1] <- exp(2 *  theta[4])
      Sigma[2, 2] <- exp(2 * theta[5])
      rho <- ilogit(Sigma[1, 2]) * 2 - 1
      Sigma[2, 1] <- Sigma[1, 2] <- rho * sqrt(Sigma[1, 1] * Sigma[2, 2])
    }
  } else {
    stop("So far we can only deal with dimension < 3")
  }
  ans <- list(alpha = alpha, Beta = Beta, Sigma = Sigma)
  ans
}

logit <- function(p) {
  log(p/(1-p))
}

ilogit <- function(p) {
  exp(p) / (1+exp(p))
}


# ------------------------------------------------------------------------
# if(FALSE) {
#   alpha <- .7
#   # p = 1
#   Beta <- matrix(1.4, 1, 1)
#   Sigma <- matrix(2.3, 1, 1)
#   # trans
#   theta <- transFunc(alpha, Beta, Sigma, TRUE)
#   theta
#   itransFunc(theta, TRUE)
#   # not trans
#   theta <- transFunc(alpha, Beta, Sigma, FALSE)
#   theta
#   itransFunc(theta, FALSE)
# 
#   # p = 2
#   Beta <- matrix(c(1.4, -.5), 1, 2)
#   Sigma <- matrix(c(1.9, -.31, -.31, 4.1), 2, 2)
#   # trans
#   theta <- transFunc(alpha, Beta, Sigma, TRUE)
#   theta
#   itransFunc(theta, TRUE)
#   # not trans
#   theta <- transFunc(alpha, Beta, Sigma, FALSE)
#   theta
#   itransFunc(theta, FALSE)
# 
# }