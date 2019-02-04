#' Fit the fractional ARMA(1,1) model.
#'
#' @template args-dX
#' @template args-dT
#' @param alpha Subdiffusion parameter, optional.
#' @template args-Tz
#' @template args-var_calc
#' @template args-dots_optim
#' @template ret-cov_vcov
#' @details The fractional ARMA(1,1) model has the form
#' \deqn{
#' \Delta X_n = \gamma \Delta X_{n-1} + \rho \Delta Z_n + (1-\theta-\rho)  \Delta Z_{n-1},
#' }
#' where \eqn{\Delta Z_n} are increments of a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by \code{trans(rho) = logit(1-rho/2)}, \code{trans(mu) = mu}, \code{\link{trans_alpha}}, and \code{\link{trans_Sigma}}.
#' @export
farma_fit <- function(dX, dT, Tz, alpha, var_calc = TRUE, ...) {
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3

  if(missing(alpha)) {
    # memory allocation
    ntheta <- 3+qq+nq
    theta_hat <- rep(NA, ntheta)
    theta_names <- c("alpha", "gamma", "rho", paste0("mu", 1:qq),
                     paste0("lambda", 1:nq))
    if(missing(Tz)) Tz <- Toeplitz(n = N)
    # profile likelihood on transformed scale
    negll.prof <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      gamma <- itrans_gamma(theta[2])
      rho <- itrans_rho(theta[3])
      Tz$setAcf(farma_acf(alpha, gamma, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.prof(suff)
    }
    # likelihood on transformed scale
    negloglik <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      gamma <- itrans_gamma(theta[2])
      rho <- itrans_rho(theta[3])
      mu <- theta[3+1:qq]
      Sigma <- itrans_Sigma(theta[3+qq+1:nq]) # default: log(D)
      Tz$setAcf(farma_acf(alpha, gamma, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    }
    # calculate MLE
    fit <- optim(fn = negll.prof, par = rep(0,3), ...)
    if(fit$convergence != 0) stop("optim did not converge.")
    theta_hat[1:3] <- fit$par # profiled parameters
    Tz$setAcf(farma_acf(itrans_alpha(theta_hat[1]), itrans_gamma(theta_hat[2]), 
                      itrans_rho(theta_hat[3]), dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    theta_hat[3+1:qq] <- suff$Beta
    theta_hat[3+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
    names(theta_hat) <- theta_names
    ans <- theta_hat # no-copy unless ans is modified
  }
  else {
    # memory allocation
    ntheta <- 2+qq+nq
    theta_hat <- rep(NA, ntheta)
    theta_names <- c("gamma", "eta", paste0("mu", 1:qq),
                     paste0("lambda", 1:nq))
    if(missing(Tz)) Tz <- Toeplitz(n = N)
    # profile likelihood on transformed scale
    negll.prof <- function(theta) {
      gamma <- itrans_gamma(theta[1])
      rho <- itrans_rho(theta[2])
      Tz$setAcf(farma_acf(alpha, gamma, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.prof(suff)
    }
    # likelihood on transformed scale
    negloglik <- function(theta) {
      gamma <- itrans_gamma(theta[1])
      rho <- itrans_rho(theta[2])
      mu <- theta[2+1:qq]
      Sigma <- itrans_Sigma(theta[2+qq+1:nq]) # default: log(D)
      Tz$setAcf(farma_acf(alpha, gamma, rho, dT, N))
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    }
    # calculate MLE
    fit <- optim(fn = negll.prof, par = rep(0,2), ...)
    if(fit$convergence != 0) stop("optim did not converge.")
    theta_hat[1:2] <- fit$par # profiled parameters
    Tz$setAcf(farma_acf(alpha, itrans_gamma(theta_hat[1]), 
                        itrans_rho(theta_hat[2]), dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    theta_hat[2+1:qq] <- suff$Beta
    theta_hat[2+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
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
