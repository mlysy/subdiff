# create farma(1,1) model object
model <- farma_model$new(Xt = Xt, dt = dt, order = c(1, 1), drift = "linear")

# evaluate loglikelihood
model$loglik(phi = c(alpha = alpha, phi = phi, rho = rho),
             mu = rep(0, ndim),
             Sigma = diag(ndim))

