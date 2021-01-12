mu <- c(1,2)
Sigma <- diag(2)
# compute the residuals
csi_resid(dX, drift = rep(dt, N) %o% mu, acf1, Sigma)
