# create fsd model object
model <- fsd_model$new(Xt = Xt, dt = dt, drift = "linear")

# evaluate loglikelihood
model$loglik(phi = c(alpha = alpha, tau = tau, sigma2 = sigma2),
             mu = rep(0, ndim),
             Sigma = diag(ndim))

