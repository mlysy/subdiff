# generates a 2-d series of fBM.
alpha <- 0.8
dt <- 1/60
N <- 1800
acf1 <- fbm_acf(alpha, dt, N)
dX <- rSnorm(2, acf = acf1)
Xt <- apply(dX, 2, cumsum)
