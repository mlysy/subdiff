context("fsd_fit")

source("subdiff-testfunctions.R")

# fsd loglikelihood
loglik_loc <- function(omega, dX, dt, Tz) {
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  alpha <- ilogit(omega[1], min = 0, max = 2)
  tau <- ilogit(omega[2], min = 0, max = 1)
  sigma2 <- exp(omega[3])
  mu <- omega[3+1:nd]
  Sigma <- itrans_Sigma(omega[3+nd+1:nq]) # default: log(D)
  acf1 <- fsd_acf(alpha, tau, sigma2, dt, N)
  Tz$set_acf(acf1)
  suff <- lmn_suff(Y = dX, X = dt, V = Tz, Vtype = "acf")
  ll <- lmn_loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  ll <- ll - log1pexp(omega[2]) - log1pexp(-omega[2]) + omega[3] # penalty term
  ll
}

ntest <- 1
test_that("MLE is at the mode of the projection plots, dynamic and localization.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(1000:2000, 1)
    dt <- runif(1)
    alpha <- runif(1, 0, 0.8)
    tau <- runif(1, 0, 1)
    sigma2 <- runif(1, 0, 0.0025)
    D <- runif(1, .9, 1.1)
    ndims <- sample(1:3, 1)
    acf1 <- D * fsd_acf(alpha, tau, sigma2, dt, N)
    dX <- as.matrix(rnormtz(n = ndims, fft = FALSE, acf = acf1))
    Xt <- apply(rbind(rnorm(ndims), dX), 2, cumsum)
    # fsd
    omega_hat <- fsd_fit(Xt, dt = dt, vcov = FALSE, ad_only = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz$new(N = N) # memory allocation
    ocheck <- optim_proj(xsol = omega_hat,
                         xnames = names(omega_hat),
                         fun = function(omega) loglik_loc(omega, dX, dt, Tz),
                         plot = FALSE, xrng = 1, npts = 50, equalize = TRUE)
    expect_lt(max_xdiff(ocheck), .05)
  }
})

## log1pe <- subdiff:::log1pe
## model <- fsd_model$new(Xt, dt)
## penalty <- function(psi) log1pe(psi[2]) + log1pe(-psi[2]) - psi[3]
## objfun1 <- model$nlp
## objfun2 <- function(psi) model$nlp(psi) + penalty(psi)

## fit1 <- optim(par = rep(0, 3), fn = objfun1)
## fit2 <- optim(par = rep(0, 3), fn = objfun2)

## optim_proj(xsol = fit1$par, fun = objfun1, maximize = FALSE)
## optim_proj(xsol = fit2$par, fun = objfun2, maximize = FALSE)

## model$get_omega(fit2$par)

# # fdy loglikelihood
# loglik1 <- function(omega, sigma2, dX, dt, Tz, penalty) {
#   N <- nrow(dX)
#   nd <- ncol(dX)
#   nq <- getq(nd)
#   alpha <- itrans_alpha(omega[1])
#   tau <- itrans_tau(omega[2])
#   mu <- omega[2+1:nd]
#   Sigma <- itrans_Sigma(omega[2+nd+1:nq]) # default: log(D)
#   acf1 <- fsd_acf(alpha, tau, sig^2, dt, N)
#   Tz$setAcf(acf1)
#   suff <- lmn.suff(Y = dX, X = dt, acf = Tz)
#   ll <- lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
#   if(penalty) {
#     ll <- ll - log1pexp(omega[2]) - log1pexp(-omega[2]) # penalty term
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
#     dt <- runif(1)
#     alpha <- runif(1, 0, 0.8)
#     tau <- runif(1, 0, 1)
#     sig <- runif(1, 0, 0.05)
#     penalty <- TRUE
#     ndims <- sample(1:2, 1)
#     acf1 <- fdyn_acf(alpha, tau, dt, N)
#     acf1[1:2] <- acf1[1:2] + sig^2 * c(2,-1)
#     dX <- as.matrix(rSnorm(n = ndims, fft = FALSE, acf = acf1))
#
#     # fdl
#     omega_hat <- fdl_fit(dX, dt, sigma2 = sig^2, vcov = FALSE, penalty = penalty,
#                          control = list(maxit = 1e6, trace = 0)) # fit MLE
#     # projection plots
#     Tz <- Toeplitz(n = N) # memory allocation
#     ocheck <- optim_proj(xsol = omega_hat,
#                          fun = function(omega) loglik1(omega, sig^2, dX, dt, Tz, penalty),
#                          plot = FALSE, xrng = .1, npts = 20, equalize = FALSE)
#     expect_lt(max.xdiff(ocheck), .05)
#   })
# })
#
# # flo likelihood
# loglik2 <- function(omega, tau, dX, dt, Tz, penalty) {
#   N <- nrow(dX)
#   nd <- ncol(dX)
#   nq <- getq(nd)
#   alpha <- itrans_alpha(omega[1])
#   sigma2 <- exp(2*omega[2])
#   mu <- omega[2+1:nd]
#   Sigma <- itrans_Sigma(omega[2+nd+1:nq]) # default: log(D)
#   acf1 <- fdyn_acf(alpha, tau, dt, N)
#   acf1[1:2] <- acf1[1:2] + sigma2 * c(2,-1)
#   Tz$setAcf(acf1)
#   suff <- lmn.suff(Y = dX, X = dt, acf = Tz)
#   ll <- lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
#   if(penalty) {
#     ll <- ll + omega[2] # penalty term
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
#     dt <- runif(1)
#     alpha <- runif(1, 0, 0.8)
#     tau <- runif(1, 0, 1)
#     sig <- runif(1, 0, 0.05)
#     penalty <- TRUE
#     ndims <- sample(1:2, 1)
#     acf1 <- fdyn_acf(alpha, tau, dt, N)
#     acf1[1:2] <- acf1[1:2] + sig^2 * c(2,-1)
#     dX <- as.matrix(rSnorm(n = ndims, fft = FALSE, acf = acf1))
#
#     # fdl
#     omega_hat <- fdl_fit(dX, dt, tau = tau, vcov = FALSE, penalty = penalty,
#                          control = list(maxit = 1e6, trace = 0)) # fit MLE
#     # projection plots
#     Tz <- Toeplitz(n = N) # memory allocation
#     ocheck <- optim_proj(xsol = omega_hat,
#                          fun = function(omega) loglik2(omega, tau, dX, dt, Tz, penalty),
#                          plot = FALSE, xrng = .1, npts = 20, equalize = FALSE)
#     expect_lt(max.xdiff(ocheck), .05)
#   })
# })
#
# # f log-likelihood
# loglik3 <- function(omega, tau, sigma2, dX, dt, Tz) {
#   N <- nrow(dX)
#   nd <- ncol(dX)
#   nq <- getq(nd)
#   alpha <- itrans_alpha(omega[1])
#   mu <- omega[1+1:nd]
#   Sigma <- itrans_Sigma(omega[1+nd+1:nq]) # default: log(D)
#   acf1 <- fdyn_acf(alpha, tau, dt, N)
#   acf1[1:2] <- acf1[1:2] + sigma2 * c(2,-1)
#   Tz$setAcf(acf1)
#   suff <- lmn.suff(Y = dX, X = dt, acf = Tz)
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
#     dt <- runif(1)
#     alpha <- runif(1, 0, 0.8)
#     tau <- runif(1, 0, 1)
#     sig <- runif(1, 0, 0.05)
#     ndims <- sample(1:2, 1)
#     acf1 <- fsd_acf(alpha, tau, sig^2, dt, N)
#     dX <- as.matrix(rSnorm(n = ndims, fft = FALSE, acf = acf1))
#
#     # fdl
#     omega_hat <- fdl_fit(dX, dt, tau = tau, sigma2 = sig^2, vcov = FALSE,
#                           control = list(maxit = 1e6, trace = 0)) # fit MLE
#     # projection plots
#     Tz <- Toeplitz(n = N) # memory allocation
#     ocheck <- optim_proj(xsol = omega_hat,
#                          fun = function(omega) loglik3(omega, tau, sig^2, dX, dt, Tz),
#                          plot = FALSE, xrng = .1, npts = 20, equalize = FALSE)
#     expect_lt(max.xdiff(ocheck), .05)
#   })
# })
