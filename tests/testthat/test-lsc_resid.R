
context("lsc_resid")

# forward transformation
lsc_fwd <- function(Z, dT, mu, acf, Sigma) {
  ed <- eigen(Sigma)
  C <- sqrt(ed$val) * t(ed$vec)
  dX <- cholZX(Z = Z, acf = acf) %*% C
  t(t(dX) + mu * dT)
}

racf <- function(N, dT, type = c("fbm", "fdl", "fma")) {
  alpha <- runif(1)
  type <- match.arg(type)
  if(type == "fbm") {
    acf1 <- fbm_acf(alpha, dT, N)
  } else if(type == "fdl") {
    tau <- runif(1)
    sig <- rexp(1)
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sig^2 * c(2,-1)
  } else if(type == "fma") {
    rho <- runif(1,-1,1)
    acf1 <- fma_acf(alpha, rho, dT, N)
  } else stop("Invalid acf type.")
  acf1
}

ntest <- 20
test_that("Z == lsc_resid(dX = lsc_fwd(Z))", {
  replicate(ntest, {
    nd <- sample(1:5, 1)
    mu <- rnorm(nd)
    Sigma <- crossprod(matrix(rnorm(nd^2),nd,nd))
    dT <- runif(1)
    N <- sample(1000:2000,1)
    Z1 <- matrix(rnorm(N*nd),N,nd)
    acf1 <- racf(N, dT, "fdl")
    dX <- lsc_fwd(Z1, dT, mu, acf1, Sigma)
    Z2 <- lsc_resid(dX, dT, mu, acf1, Sigma)
    expect_equal(max(abs(Z1-Z2)), 0)
  })
})
