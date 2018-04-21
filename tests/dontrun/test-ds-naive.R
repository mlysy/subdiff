require(subdiff)
source("D:/GitHub/SubDiff/tests/dontrun/mle-check.R")

# paremeters
alpha <- 0.8
N <- 1000
dT <- 1/60
Beta <- matrix(c(1.4, 2.6), 1, 2)
Sigma <- matrix(c(2, 0, 0, 3), 2, 2)

# generating time series
acf1 <- fbm.acf(alpha, dT, N)
dY <- rbind(0, rSnorm(n = 2, acf = acf1))
Yt <- apply(dY, 2, cumsum)
tSeq <- 0:N + 1
Xt <- as.matrix(tSeq * dT)
Yt <- as.matrix(Yt) %*% Sigma + Xt %*% Beta

plot(Yt[, 1], type = "l")

#  ------------------------------------------------------------------------


ds <- 5
N2 <- floor((N+1)/ds) - 1
Tz2 <- Toeplitz(n = N2)

theta <- downsample.naive(Xt = Yt, dT = dT, Tz = Tz2, ds = ds)

comp.test3 <- function(theta) {
  suff <- itransFunc(theta, FALSE)
  Tz2$setAcf(fbm.acf(suff$alpha, ds*dT, N2))
  ans <- composite.full(Y = downSample(Yt, ds), X = ds*dT, suff$Beta, suff$Sigma, acf = Tz2, ds = 1)
  ans
}

mle.check(comp.test3, theta)
