
context("far_fit")

source("fit-functions.R")

# fma profile likelihood
profll <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- itrans_alpha(theta[1])
  rho <- itrans_rho(theta[2])
  mu <- theta[2+1:nd]
  Sigma <- itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
  Tz$setAcf(fbm_acf(alpha, dT, N))
  dY <- ar_resid(dX, rho, denom = FALSE)
  suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
  ll <- lmn.prof(suff) # - log(1 - rho) * N * nd
  # penalty
  #ll - log1pexp(theta[2]) - log1pexp(-theta[2])
  ll
}


# fma loglikelihood
loglik <- function(theta, dX, dT, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- itrans_alpha(theta[1])
  rho <- itrans_rho(theta[2])
  mu <- (1-rho) * theta[2+1:nd]
  Sigma <- (1-rho)^2 * itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
  Tz$setAcf(fbm_acf(alpha, dT, N))
  dY <- ar_resid(dX, rho, denom = FALSE)
  suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
  lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff) # - log(1 - rho) * N * nd
}

ntest <- 10
test_that("MLE is at the mode of the projection plots.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  replicate(n = ntest, {
    # simulate data
    N <- sample(1000:2000, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 2)
    rho <- runif(1, 0, 1)
    ndims <- sample(1:2, 1)
    dX <- as.matrix(rSnorm(n = ndims,
                           acf = fbm_acf(alpha, dT, N)))
    dY <- apply((1-rho)*dX, 2, function(x) {
      as.numeric(filter(x, filter = rho, method = "recursive"))
    })
    theta_hat <- far_fit(dY, dT, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat[1:2],
                         fun = function(theta) profll(theta, dY, dT, Tz),
                         plot = FALSE,
                         xrng = .1,
                         #xrng = cbind(c(-5,5), c(-10, 30)),
                         npts = 20, equalize = FALSE)
    expect_lt(max.xdiff(ocheck), .05)
  })
})


##     # simulate data
##     N <- sample(1000:2000, 1)
##     dT <- runif(1)
##     alpha <- runif(1, 0, 2)
##     rho <- runif(1, 0, 1)
##     ndims <- sample(1:2, 1)
##     dX <- as.matrix(rSnorm(n = ndims,
##                            acf = fbm_acf(alpha, dT, N)))
##     dY <- apply((1-rho)*dX, 2, function(x) {
##       as.numeric(filter(x, filter = rho, method = "recursive"))
##     })
## theta_hat <- far_fit(dY, dT, var_calc = TRUE) # fit MLE
## theta_hat2 <- far_fit2(dY, dT, var_calc = TRUE)
##     # projection plots
##     Tz <- Toeplitz(n = N) # memory allocation
##     ocheck <- optim_proj(xsol = theta_hat[1:2],
##                          fun = function(theta) profll(theta, dY, dT, Tz),
##                          plot = TRUE,
##                          #xrng = 1,
##                          xrng = cbind(c(-5,5), c(-10, 10)),
##                          npts = 100, equalize = TRUE)
