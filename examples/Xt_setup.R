# simulate a 2D fbm trajectory.
alpha <- 0.8
dt <- 1/60
N <- 1800
acf <- fbm_acf(alpha, dt, N)
dX <- SuperGauss::rnormtz(n = 2, acf = acf)
Xt <- apply(dX, 2, cumsum)

