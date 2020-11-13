
context("csi_resid")

# forward transformation
csi_fwd <- function(Z, dt, mu, acf, Sigma) {
  ed <- eigen(Sigma)
  C <- sqrt(ed$val) * t(ed$vec)
  dX <- SuperGauss::cholZX(Z = Z, acf = acf) %*% C
  t(t(dX) + mu * dt)
}

racf <- function(N, dt, type = c("fbm", "fsd", "farma")) {
  alpha <- runif(1)
  type <- match.arg(type)
  if(type == "fbm") {
    acf1 <- fbm_acf(alpha, dt, N)
  } else if(type == "fsd") {
    tau <- runif(1)
    sig <- rexp(1)
    acf1 <- fsd_acf(alpha, tau, sig^2, dt, N)
  } else if(type == "farma") {
    phi <- rho <- runif(1,-1,1)
    acf1 <- farma_acf(alpha, phi, rho, dt, N)
  } else stop("Invalid acf type.")
  acf1
}

ntest <- 20
test_that("Z == csi_resid(dX = csi_fwd(Z))", {
  replicate(ntest, {
    nd <- sample(1:5, 1)
    mu <- rnorm(nd)
    Sigma <- crossprod(matrix(rnorm(nd^2),nd,nd))
    dt <- runif(1)
    N <- sample(1000:2000,1)
    Z1 <- matrix(rnorm(N*nd),N,nd)
    acf1 <- racf(N, dt, "fsd")
    dX <- csi_fwd(Z1, dt, mu, acf1, Sigma)
    Z2 <- csi_resid(dX, dt, mu, acf1, Sigma)
    expect_equal(max(abs(Z1-Z2)), 0)
  })
})
