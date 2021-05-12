# simulate data from the fbm model
alpha <- .8

dt <- 1/60
N <- 1800
ndim <- 2

Xt <- csi_sim(drift = matrix(0, N-1, ndim),
              acf = fbm_acf(alpha, dt, N-1),
              Sigma = diag(ndim),
              X0 = rep(0, ndim))


