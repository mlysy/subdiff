require(subdiff)
require(mniw)
source("D:/GitHub/SubDiff/tests/dontrun/mle-check.R")

# paremeters
alpha <- 0.8
N <- 1800
dT <- 1/60
Beta <- matrix(c(1.4, 2.6), 1, 2)
Sigma <- matrix(c(2, 0, 0, 3), 2, 2)
ds <- 5

# comparison function
test_full_loglik <- function(alpha, mu, Sigma, Xt, dT, ds) {
  ll <- 0
  for(ii in 1:ds) {
    X <- downSample(Yt = Xt, ds = ds, pos = ii)
    X <- apply(X, 2, diff)
    n <- nrow(X)
    Mu <- rep(dT*ds, n) %*% mu
    V <- toeplitz(fbm.acf(alpha, dT*ds, n))
    ll <- ll + dMNorm(X, Mu = Mu, RowV = V, ColV = Sigma, log = TRUE)
  }
  ll
}

for(ii in 1:3) {
  Yt <- matrix(rnorm(2*N), N, 2)
  acf1 <- fbm.acf(alpha = alpha, dT = ds*dT, N = floor(N/ds)-1)
  Tz <- Toeplitz(acf = acf1)
  ll1 <- composite.full(Y = Yt, X = dT, Beta = Beta, Sigma = Sigma, acf = Tz, ds = ds)
  ll2 <- test_full_loglik(alpha = alpha, mu = Beta, Sigma = Sigma, Xt = Yt, dT = dT, ds = ds)
  print(ll1 - ll2)
}


