Tz_fit <- function(ntheta, acf_func, dX, dT, Tz, var_calc, penalty = FALSE) {
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  theta_hat <- rep(NA, ntheta+qq+nq)
  
  # profile likelihood on transformed scale
  negll.prof <- function(theta) {
    Tz$setAcf(acf_func(theta))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    nlp <- -lmn.prof(suff)
    if(penalty) nlp <- nlp + log1pe(theta[2]) + log1pe(-theta[2]) - theta[3]
    nlp
  }

  # calculate MLE
  if(ntheta == 1) {
    theta_hat[1] <- optimize(f = negll.prof, interval = c(-5, 5))$minimum
  } else {
    fit <- optim(fn = negll.prof, par = rep(0,ntheta))
    if(fit$convergence != 0) warning("optim did not converge.")
    theta_hat[1:ntheta] <- fit$par
  }
  
  Tz$setAcf(acf_func(theta_hat[1:ntheta]))
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[ntheta+1:qq] <- suff$Beta
  theta_hat[ntheta+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
  ans <- theta_hat
  
  if(var_calc) {
    # likelihood on transformed scale
    negloglik <- function(theta) {
      mu <- theta[ntheta+1:qq]
      Sigma <- itrans_Sigma(theta[ntheta+qq+1:nq])
      Tz$setAcf(acf_func(theta[1:ntheta]))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      nlp <- -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
      if(penalty) nlp <- nlp + log1pe(theta[2]) + log1pe(-theta[2]) - theta[3]
      nlp
    }
    # variance estimate
    V_hat <- hessian(negloglik, x = theta_hat)
    V_hat <- solveV(V_hat)
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}

