
context("fdl_fit")

source("fit-functions.R")

# fdl loglikelihood
loglik <- function(theta, dX, dT, Tz, penalty) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- itrans_alpha(theta[1])
  tau <- itrans_tau(theta[2])
  sigma2 <- exp(2*theta[3])
  mu <- theta[3+1:nd]
  Sigma <- itrans_Sigma(theta[3+nd+1:nq]) # default: log(D)
  acf1 <- fdyn_acf(alpha, tau, dT, N)
  acf1[1:2] <- acf1[1:2] + sigma2 * c(2,-1)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  ll <- lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  ll <- ll - log1pexp(theta[2]) - log1pexp(-theta[2]) + theta[3] # penalty term
  ll
}

ntest <- 10
test_that("MLE is at the mode of the projection plots, dynamic and localization.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  replicate(n = ntest, {
    # simulate data
    N <- sample(1000:2000, 1)
    dT <- runif(1)
    alpha <- runif(1, 0, 0.8)
    tau <- runif(1, 0, 1)
    sig <- runif(1, 0, 0.05)
    penalty <- TRUE
    ndims <- sample(1:2, 1)
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sig^2 * c(2,-1)
    dX <- as.matrix(rSnorm(n = ndims, fft = FALSE, acf = acf1))
    
    # fdl
    theta_hat <- floc_fit(dX, dT, var_calc = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz(n = N) # memory allocation
    ocheck <- optim_proj(xsol = theta_hat,
                         fun = function(theta) loglik(theta, dX, dT, Tz, penalty),
                         plot = FALSE, xrng = .1, npts = 20, equalize = FALSE)
    
    expect_lt(max.xdiff(ocheck), .05)
  })
})

# # fdy loglikelihood
# loglik1 <- function(theta, sigma2, dX, dT, Tz, penalty) {
#   N <- nrow(dX)
#   nd <- ncol(dX)
#   nq <- getq(nd)
#   alpha <- itrans_alpha(theta[1])
#   tau <- itrans_tau(theta[2])
#   mu <- theta[2+1:nd]
#   Sigma <- itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
#   acf1 <- fdyn_acf(alpha, tau, dT, N)
#   acf1[1:2] <- acf1[1:2] + sigma2 * c(2,-1)
#   Tz$setAcf(acf1)
#   suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
#   ll <- lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
#   if(penalty) {
#     ll <- ll - log1pexp(theta[2]) - log1pexp(-theta[2]) # penalty term
#   }
#   ll
# }
# 
# ntest <- 10
# test_that("MLE is at the mode of the projection plots, dynamic.", {
#   skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
#               "Package \"optimCheck\" required to run this test.")
#   require(optimCheck)
#   replicate(n = ntest, {
#     # simulate data
#     N <- sample(1000:2000, 1)
#     dT <- runif(1)
#     alpha <- runif(1, 0, 0.8)
#     tau <- runif(1, 0, 1)
#     sig <- runif(1, 0, 0.05)
#     penalty <- TRUE
#     ndims <- sample(1:2, 1)
#     acf1 <- fdyn_acf(alpha, tau, dT, N)
#     acf1[1:2] <- acf1[1:2] + sig^2 * c(2,-1)
#     dX <- as.matrix(rSnorm(n = ndims, fft = FALSE, acf = acf1))
#     
#     # fdl
#     theta_hat <- fdl_fit(dX, dT, sigma2 = sig^2, var_calc = FALSE, penalty = penalty,
#                          control = list(maxit = 1e6, trace = 0)) # fit MLE
#     # projection plots
#     Tz <- Toeplitz(n = N) # memory allocation
#     ocheck <- optim_proj(xsol = theta_hat,
#                          fun = function(theta) loglik1(theta, sig^2, dX, dT, Tz, penalty),
#                          plot = FALSE, xrng = .1, npts = 20, equalize = FALSE)
#     expect_lt(max.xdiff(ocheck), .05)
#   })
# })
# 
# # flo likelihood
# loglik2 <- function(theta, tau, dX, dT, Tz, penalty) {
#   N <- nrow(dX)
#   nd <- ncol(dX)
#   nq <- getq(nd)
#   alpha <- itrans_alpha(theta[1])
#   sigma2 <- exp(2*theta[2])
#   mu <- theta[2+1:nd]
#   Sigma <- itrans_Sigma(theta[2+nd+1:nq]) # default: log(D)
#   acf1 <- fdyn_acf(alpha, tau, dT, N)
#   acf1[1:2] <- acf1[1:2] + sigma2 * c(2,-1)
#   Tz$setAcf(acf1)
#   suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
#   ll <- lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
#   if(penalty) {
#     ll <- ll + theta[2] # penalty term
#   }
#   ll
# }
# 
# ntest <- 10
# test_that("MLE is at the mode of the projection plots, localization.", {
#   skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
#               "Package \"optimCheck\" required to run this test.")
#   require(optimCheck)
#   replicate(n = ntest, {
#     # simulate data
#     N <- sample(1000:2000, 1)
#     dT <- runif(1)
#     alpha <- runif(1, 0, 0.8)
#     tau <- runif(1, 0, 1)
#     sig <- runif(1, 0, 0.05)
#     penalty <- TRUE
#     ndims <- sample(1:2, 1)
#     acf1 <- fdyn_acf(alpha, tau, dT, N)
#     acf1[1:2] <- acf1[1:2] + sig^2 * c(2,-1)
#     dX <- as.matrix(rSnorm(n = ndims, fft = FALSE, acf = acf1))
#     
#     # fdl
#     theta_hat <- fdl_fit(dX, dT, tau = tau, var_calc = FALSE, penalty = penalty,
#                          control = list(maxit = 1e6, trace = 0)) # fit MLE
#     # projection plots
#     Tz <- Toeplitz(n = N) # memory allocation
#     ocheck <- optim_proj(xsol = theta_hat,
#                          fun = function(theta) loglik2(theta, tau, dX, dT, Tz, penalty),
#                          plot = FALSE, xrng = .1, npts = 20, equalize = FALSE)
#     expect_lt(max.xdiff(ocheck), .05)
#   })
# })
# 
# # f log-likelihood
# loglik3 <- function(theta, tau, sigma2, dX, dT, Tz) {
#   N <- nrow(dX)
#   nd <- ncol(dX)
#   nq <- getq(nd)
#   alpha <- itrans_alpha(theta[1])
#   mu <- theta[1+1:nd]
#   Sigma <- itrans_Sigma(theta[1+nd+1:nq]) # default: log(D)
#   acf1 <- fdyn_acf(alpha, tau, dT, N)
#   acf1[1:2] <- acf1[1:2] + sigma2 * c(2,-1)
#   Tz$setAcf(acf1)
#   suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
#   ll <- lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
#   ll
# }
# 
# ntest <- 10
# test_that("MLE is at the mode of the projection plots.", {
#   skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
#               "Package \"optimCheck\" required to run this test.")
#   require(optimCheck)
#   replicate(n = ntest, {
#     # simulate data
#     N <- sample(1000:2000, 1)
#     dT <- runif(1)
#     alpha <- runif(1, 0, 0.8)
#     tau <- runif(1, 0, 1)
#     sig <- runif(1, 0, 0.05)
#     ndims <- sample(1:2, 1)
#     acf1 <- fdyn_acf(alpha, tau, dT, N)
#     acf1[1:2] <- acf1[1:2] + sig^2 * c(2,-1)
#     dX <- as.matrix(rSnorm(n = ndims, fft = FALSE, acf = acf1))
#     
#     # fdl
#     theta_hat <- fdl_fit(dX, dT, tau = tau, sigma2 = sig^2, var_calc = FALSE,
#                           control = list(maxit = 1e6, trace = 0)) # fit MLE
#     # projection plots
#     Tz <- Toeplitz(n = N) # memory allocation
#     ocheck <- optim_proj(xsol = theta_hat,
#                          fun = function(theta) loglik3(theta, tau, sig^2, dX, dT, Tz),
#                          plot = FALSE, xrng = .1, npts = 20, equalize = FALSE)
#     expect_lt(max.xdiff(ocheck), .05)
#   })
# })
