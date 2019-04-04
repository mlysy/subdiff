# generating time series
dT <- 1/60
N <- 1800
acf1 <- fbm_acf(0.8, dT, N)
dX <- rSnorm(n = 2, acf = acf1)
Tz <- Toeplitz(n = N)