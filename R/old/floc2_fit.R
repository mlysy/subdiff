#' Fit the fBM model with dynamic and static errors.
#'
#' @name floc2_fit
#' @template args-dX
#' @template args-dT
#' @param tau Scalar between 0 and 1, indicating the percentage of time for which the camera shutter is open in the dynamic error model.  Estimated if missing. See Details.
#' @param sigma2 Magnitude of static error.  Estimated if missing. See Details.
#' @template args-Tz
#' @template args-var_calc
#' @param penalty logic, employ a small penalty on `(tau, sigma2)` when `TRUE`
#' @template ret-cov_vcov 
#' @details The fBM + dynamic and localization error (fdl) model has the form
#' \deqn{
#' X_n = 1/\tau \int_0^\tau Z_(n+s) ds + \sigma e_{n},
#' }
#' where \eqn{Z_n} is a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by `trans(tau) = logit(tau)`, `trans(sigma) = log(sigma)`, `trans(mu) = mu`, [trans_alpha()], and [trans_Sigma()].
#' When put `tau` = 0, this model becomes fractional localization error model.
#' When put `sigma2` = 0, this model becomes fractional dynamic error model.
#' When put `(tau, sigma2)` = 0, this model becomes fractional Brownian model.
#' 
#' @export
floc2_fit <- function(dX, dT, tau, Tz, var_calc = TRUE) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  
  theta_names <- c("alpha", "sigma2", paste0("mu", 1:qq),
                   paste0("lambda", 1:nq))
  # acf function on transformed scale
  acf_func <- function(theta) {
    floc_acf(itrans_alpha(theta[1]), tau, exp(2*theta[2]), dT, N)
  }
  
  # estimation
  ans <- Tz2_fit(2, acf_func, dX, dT, Tz, var_calc, penalty = TRUE)
  
  # add names
  if(var_calc) {
    names(ans$coef) <- colnames(ans$vcov) <- 
      rownames(ans$vcov) <- theta_names
  } else {
    names(ans) <- theta_names
  }
  ans
}


Tz2_fit <- function(ntheta, acf_func, dX, dT, Tz, var_calc, penalty = FALSE) {
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  theta_hat <- rep(NA, ntheta+qq+nq)
  
  # profile likelihood on transformed scale
  negll.prof <- function(theta) {
    Tz$setAcf(acf_func(theta))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    nlp <- -lmn.prof(suff)
    if(penalty) nlp <- nlp - theta[2]
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
      if(penalty) nlp <- nlp - theta[2]
      nlp
    }
    # variance estimate
    V_hat <- hessian(negloglik, x = theta_hat)
    V_hat <- solveV(V_hat)
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}

