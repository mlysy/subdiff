# simulate data from the fsd model
alpha <- .8
tau <- .1
sigma2 <- .01

dt <- 1/60
N <- 1800
ndim <- 2

Xt <- csi_sim(drift = matrix(0, N-1, ndim),
              acf = fsd_acf(alpha, tau, sigma2, dt, N-1),
              Sigma = diag(ndim),
              X0 = rep(0, ndim))


