
context("fbm_fit")

source("fit-functions.R")

# fbm loglikelihood
loglik <- function(omega, dX, dt, Tz,
                   drift = c("linear", "none", "quadratic")) {
  drift <- match.arg(drift)
  N <- nrow(dX)
  nd <- ncol(dX)
  nq <- getq(nd)
  dr <- switch(drift,
               none = drift_none(0, dt, N),
               linear = drift_linear(0, dt, N),
               quadratic = drift_quadratic(0, dt, N))
  mobj <- fbm_model$new(dX[1,,drop=FALSE], 1, drift = drift)
  theta <- mobj$itrans_full(omega)
  ## alpha <- ilogit(theta[1], min = 0, max = 2)
  ## mu <- theta[1+1:nd]
  ## Sigma <- itrans_Sigma(theta[1+nd+1:nq]) # default: log(D)
  Tz$set_acf(fbm_acf(theta$phi, dt, N))
  suff <- lmn_suff(Y = dX, X = dr, V = Tz, Vtype = "acf")
  lmn_loglik(Beta = theta$mu, Sigma = theta$Sigma, suff = suff)
}

test_that("MLE is at the mode of the projection plots.", {
  skip_if_not(requireNamespace("optimCheck", quietly = TRUE),
              "Package \"optimCheck\" required to run this test.")
  require(optimCheck)
  test_cases <- expand.grid(drift = c("linear", "none", "quadratic"),
                            ndims = 1:2)
  ntest <- nrow(test_cases)
  for(ii in 1:ntest) {
    # simulate data
    N <- sample(500:1000, 1)
    dt <- runif(1)
    alpha <- runif(1, 0, 2)
    D <- runif(1, .9, 1.1)
    ndims <- test_cases$ndims[ii]
    dX <- as.matrix(rnormtz(n = ndims,
                            acf = D * fbm_acf(alpha, dt, N)))
    drift <- as.character(test_cases$drift[ii])
    omega_hat <- fbm_fit(dX, dt, drift = drift,
                         vcov = FALSE, ad_only = FALSE) # fit MLE
    # projection plots
    Tz <- Toeplitz$new(N = N) # memory allocation
    ocheck <- optim_proj(xsol = omega_hat,
                         fun = function(omega) loglik(omega, dX, dt, Tz, drift),
                         plot = FALSE, xrng = .05, npts = 20)
    expect_lt(max_xdiff(ocheck), .01)
  }
})
