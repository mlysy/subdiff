require(subdiff)
source("D:/GitHub/SubDiff/tests/dontrun/mle-check.R")
source("D:/GitHub/SubDiff/tests/dontrun/composite-full.R")

alpha <- 0.8
N <- 1800
dt <- 1/60
acf1 <- fbm.acf(alpha, dt, N)
dY <- rSnorm(n = 1, acf = acf1)
Yt <- c(0, cumsum(dY))
tSeq <- 0:N + 1
Xt <- as.matrix(tSeq)
Yt <- as.matrix(Xt * 1.4 + Yt)
Toep1 <- Toeplitz(acf = acf1)

suff1 <- composite.suff(Y = Yt, X = Xt, acf = Toep1, ds = 1)
suff1$Betahat
suff1$S / N

N2 <- floor((N+1)/2) - 1
acf2 <- fbm.acf(alpha, dt * 2, N2)
Toep2 <- Toeplitz(acf = acf2)

suff2 <- composite.suff(Y = Yt, X = Xt, acf = Toep2, ds = 2)
suff2$Betahat
suff2$S / N2 / 2

N3 <- floor((N+1)/3) - 1
acf3 <- fbm.acf(alpha, dt * 3, N3)
Toep3 <- Toeplitz(acf = acf3)

suff3 <- composite.suff(Y = Yt, X = Xt, acf = Toep3, ds = 3)
suff3$Betahat
suff3$S / N3 / 3

N4 <- floor((N+1)/4) - 1
acf4 <- fbm.acf(alpha, dt * 4, N4)
Toep4 <- Toeplitz(acf = acf4)

suff4 <- composite.suff(Y = Yt, X = Xt, acf = Toep4, ds = 4)
suff4$Betahat
suff4$S / N4 / 4

comp.test <- function(theta) {
  Beta <- as.matrix(theta[1])
  Sigma <- as.matrix(theta[2])
  ans <- composite.full(Y = Yt, X = Xt, Beta, Sigma, acf = Toep4, ds = 4)
  ans
}

mle.check(comp.test, c(as.numeric(suff4$Betahat), as.numeric(suff4$S/N4/4)))
