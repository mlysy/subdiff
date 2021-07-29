# compute the residuals
Z <- csi_resid(dX = apply(Xt, 2, diff),
               drift = matrix(0, N-1, ndim),
               acf = fbm_acf(alpha, dt, N-1),
               Sigma = diag(ndim))
