# ok how to efficiently work with subdiff models.

# characteristics:
# - all 2D location-scale models
# - each have parameters alpha and D

# generic methods:
# - fit: does model fitting
# - coef & vcov: extracts (named) coefficients and estimated variance matrix.
#   these should potentially return:
#   * alpha, D
#   * alpha, log(D)
#   * alpha, mu, Sigma
#   * alpha, mu, LSR
#   * logit(alpha), ...?
# - residuals: for cholesky time decomp, eigen space decomp

# order of variables: alpha, ..., mu, Sigma.

# ok parametrization:
# Sigma <-> LSR: nu = log(D), kappa = log(Sig1/Sig2), omega = logit(.5*rho + .5)
# alpha <-> gamma = logit(alpha/2)

# so for example:
# coef(M, type = c("aD", "full"), trans = FALSE, itrans_fun)

# aD: only (alpha, D), or (gamma, nu) if trans = TRUE
# full: (alpha, ..., mu, sig1, sig2, rho) or (gamma, ..., mu, nu, kappa, omega)
# itrans_fun specifies the inverse transformation to be used

# vcov(M, type = c("aD", "full"), trans = FALSE, itrans_fun)

# ok transform functions.  only for q = 1 and q = 2



