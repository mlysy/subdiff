# generating 2-d stationary time series
alpha <- 0.8
dT <- 1/60
N <- 1800
acf1 <- fbm_acf(alpha, dT, N)
dX <- rSnorm(2, acf = acf1)
# Toeplitz object for internal computation
Tz <- Toeplitz(n = N)