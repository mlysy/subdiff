# create fbm model object
model <- fbm_model$new(Xt = Xt, dt = dt, drift = "linear")

# evaluate loglikelihood
model$loglik(phi = c(alpha = alpha),
             mu = rep(0, ndim),
             Sigma = diag(ndim))

