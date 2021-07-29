# simulate data from a farma(1,1) model
alpha <- .8
phi <- .1
rho <- .1

dt <- 1/60
N <- 1800
ndim <- 2

Xt <- csi_sim(drift = matrix(0, N-1, ndim),
              acf = farma_acf(alpha, phi = phi, rho = rho, dt, N-1),
              Sigma = diag(ndim),
              X0 = rep(0, ndim))


