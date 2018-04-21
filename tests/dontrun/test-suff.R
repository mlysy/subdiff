require(subdiff)
source("D:/GitHub/SubDiff/tests/dontrun/mle-check.R")

# paremeters
alpha <- 0.8
N <- 100
dT <- 1/60
Beta <- matrix(c(1.4, 2.6), 1, 2)
Sigma <- matrix(c(2, 0, 0, 3), 2, 2)

# generating time series
acf1 <- fbm_acf(alpha, dT, N)
dY <- rbind(0, rSnorm(n = 2, acf = acf1))
Yt <- apply(dY, 2, cumsum)
tSeq <- 0:N + 1
Xt <- as.matrix(tSeq * dT)
Yt <- as.matrix(Yt) %*% Sigma + Xt %*% Beta

plot(Yt[, 1], type = "l")

#  ------------------------------------------------------------------------

ds <- 4
N2 <- floor((N+1)/ds) - 1
Tz2 <- Toeplitz(n = N2)

alpha2 <- optimize(f = test.fbm.prof.cl, interval = c(0, 2), Xt = Yt, dT = dT, ds = ds, Tz = Tz2, 
                  maximum = TRUE)$maximum
acf2 <- fbm_acf(alpha2, ds*dT, N2)
Tz2$setAcf(acf2)
suff2 <- composite.suff(Y = Yt, X = dT, acf = Tz2, ds = ds)
suff2$S <- suff2$S/N2/ds
theta <- c(alpha2, suff2$Beta[1], suff2$Beta[2], suff2$S[1, 1], suff2$S[2, 2], suff2$S[2, 1])

comp.test2 <- function(theta) {
  Beta <- matrix(c(theta[2], theta[3]), 1, 2)
  Sigma <- matrix(c(theta[4], theta[6], theta[6], theta[5]), 2, 2)
  acf2 <- fbm_acf(theta[1], ds*dT, N2)
  Tz2$setAcf(acf2)
  ans <- composite.full(Y = Yt, X = dT, Beta, Sigma, acf = Tz2, ds = ds)
  ans
}

mle.check(comp.test2, theta)
