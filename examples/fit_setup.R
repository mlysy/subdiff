# generating 2-d stationary time series
alpha <- 0.8
dt <- 1/60
N <- 1800
acf1 <- fbm_acf(alpha, dt, N)
dX <- rSnorm(2, acf = acf1)
# Toeplitz object for internal computation
Tz <- Toeplitz(n = N)
