# generates a 2-d series of fBM.
alpha <- 0.8
dT <- 1/60
N <- 1800
acf1 <- fbm_acf(alpha, dT, N)
dX <- rSnorm(2, acf = acf1)
Xt <- apply(dX, 2, cumsum)
