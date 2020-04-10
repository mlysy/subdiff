#--- Implementation of subdiff model_object and model_fit ----------------------

# essentially memory allocation.
fmod <- fbm_model(dX, drift = "linear")

# usage:

# basics
fmod$acf(phi) # original parametrization
fmod$drift(phi)
fmod$mu_hat(phi)
fmod$Sigma_hat(phi)
# parameter conversions
fmod$trans(mu, Sigma, phi)
fmod$itrans(psi)

# model fitting
fmod$nlp(psi) # computational parametrization
fmod$nlp_grad(psi) # if applicable
fmod$loglik(psi)

# model residuals
fmod$resid(mu, Sigma, phi) # original scale???

# change data? maybe not worth it...
fmod$set_dX(dX)

