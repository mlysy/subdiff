
context("csi_resid")

# forward transformation
csi_fwd <- function(Z, dT, mu, acf, Sigma) {
  ed <- eigen(Sigma)
  C <- sqrt(ed$val) * t(ed$vec)
  dX <- SuperGauss::cholZX(Z = Z, acf = acf) %*% C
  t(t(dX) + mu * dT)
}

racf <- function(N, dT, type = c("fbm", "floc", "farma")) {
  alpha <- runif(1)
  type <- match.arg(type)
  if(type == "fbm") {
    acf1 <- fbm_acf(alpha, dT, N)
  } else if(type == "floc") {
    tau <- runif(1)
    sig <- rexp(1)
    acf1 <- floc_acf(alpha, tau, sig^2, dT, N)
  } else if(type == "farma") {
    phi <- rho <- runif(1,-1,1)
    acf1 <- farma_acf(alpha, phi, rho, dT, N)
  } else stop("Invalid acf type.")
  acf1
}

ntest <- 20
test_that("Z == csi_resid(dX = csi_fwd(Z))", {
  replicate(ntest, {
    nd <- sample(1:5, 1)
    mu <- rnorm(nd)
    Sigma <- crossprod(matrix(rnorm(nd^2),nd,nd))
    dT <- runif(1)
    N <- sample(1000:2000,1)
    Z1 <- matrix(rnorm(N*nd),N,nd)
    acf1 <- racf(N, dT, "floc")
    dX <- csi_fwd(Z1, dT, mu, acf1, Sigma)
    Z2 <- csi_resid(dX, dT, mu, acf1, Sigma)
    expect_equal(max(abs(Z1-Z2)), 0)
  })
})
