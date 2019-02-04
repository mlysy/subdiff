#' Fit the fractional Moving-average model.
#'
#' @template args-dX
#' @template args-dT
#' @param nlag Number of lags (see Details).
#' @param alpha Subdiffusion parameter, optional.
#' @template args-Tz
#' @template args-var_calc
#' @template args-dots_optim
#' @template ret-cov_vcov
#' @details The fractional Moving-average model has the form
#' \deqn{
#' \Delta X_n = (1-\sum_{i=1}^p \rho_i) \Delta Z_n + \sum_{i=1}^p \rho_{i} \Delta Z_{n-i},
#' }
#' where \eqn{\Delta Z_n} are increments of a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by \code{trans(rho) = logit(1-rho/2)}, \code{trans(mu) = mu}, \code{\link{trans_alpha}}, and \code{\link{trans_Sigma}}.
#' @export
fma_fit <- function(dX, dT, nlag, Tz, alpha, var_calc = TRUE, ...) {
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  if(missing(nlag)) nlag <- 1
  
  if(missing(alpha)) {
    # memory allocation
    ntheta <- 1+nlag+qq+nq
    theta_hat <- rep(NA, ntheta)
    theta_names <- c("gamma", paste0("rho", 1:nlag), paste0("mu", 1:qq),
                     paste0("lambda", 1:nq))
    if(missing(Tz)) Tz <- Toeplitz(n = N)
    # profile likelihood on transformed scale
    negll.prof <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      rho <- theta[1+1:nlag]
      Tz$setAcf(fma_acf(alpha, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.prof(suff)
    }
    # likelihood on transformed scale
    negloglik <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      rho <- theta[1+1:nlag]
      mu <- theta[1+nlag+1:qq]
      Sigma <- itrans_Sigma(theta[1+nlag+qq+1:nq]) # default: log(D)
      Tz$setAcf(fma_acf(alpha, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    }
    # calculate MLE
    fit <- optim(fn = negll.prof, par = rep(0,nlag+1), ...)
    if(fit$convergence != 0) stop("optim did not converge.")
    theta_hat[1+0:nlag] <- fit$par # profiled parameters
    Tz$setAcf(fma_acf(itrans_alpha(theta_hat[1]), theta_hat[1+1:nlag], dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    theta_hat[1+nlag+1:qq] <- suff$Beta
    theta_hat[1+nlag+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
    names(theta_hat) <- theta_names
    ans <- theta_hat # no-copy unless ans is modified
  }
  else {
    # memory allocation
    ntheta <- nlag+qq+nq
    theta_hat <- rep(NA, ntheta)
    theta_names <- c(paste0("eta", 1:nlag), paste0("mu", 1:qq),
                     paste0("lambda", 1:nq))
    if(missing(Tz)) Tz <- Toeplitz(n = N)
    # profile likelihood on transformed scale
    negll.prof <- function(theta) {
      rho <- theta[1:nlag]
      Tz$setAcf(fma_acf(alpha, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.prof(suff)
    }
    # likelihood on transformed scale
    negloglik <- function(theta) {
      rho <- theta[1:nlag]
      mu <- theta[nlag+1:qq]
      Sigma <- itrans_Sigma(theta[nlag+qq+1:nq]) # default: log(D)
      Tz$setAcf(fma_acf(alpha, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    }
    # calculate MLE
    fit <- optimize(f = negll.prof, interval = c(-10, 10))
    theta_hat[0:nlag] <- fit$minimum # profiled parameters
    Tz$setAcf(fma_acf(alpha, theta_hat[1:nlag], dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    theta_hat[nlag+1:qq] <- suff$Beta
    theta_hat[nlag+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
    names(theta_hat) <- theta_names
    ans <- theta_hat # no-copy unless ans is modified
  }
  if(var_calc) {
    # variance estimate
    V_hat <- hessian(negloglik, x = theta_hat)
    V_hat <- solveV(V_hat)
    colnames(V_hat) <- theta_names
    rownames(V_hat) <- theta_names
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}

