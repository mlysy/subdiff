#' Fit the fractional dynamic localization error model.
#'
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param Tz Optional Toeplitz matrix for intermediate calculations.
#' @param type One of "dynamic localization", "dynamic" or "localization".
#' @param var_calc If \code{TRUE}, also estimate variance matrix.
#' @param ... Additional \code{control} arguments to \code{stats::optim}.
#' @return Vector of coefficients and possibly variance matrix on the transformed scale (see Details).
#' @details The fractional AR(1) model has the form
#' \deqn{
#' X_n = 1/\tau \int_0^\tau Z_(n+s) ds + \sigma e_{n},
#' }
#' where \eqn{Z_n} is a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed 
#' scale defined by \code{trans(rho) = logit(1-rho/2)}, \code{trans(mu) = mu}, \code{\link{trans_alpha}}, and 
#' \code{\link{trans_Sigma}}.
#' @export
fdl_fit <- function(dX, dT, Tz, type = "dynamic localization", var_calc = TRUE, ...) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  nq <- if(q == 1) 1 else 3
  if(type == "dynamic localization") {
    ntheta <- 3+q+nq
    theta_hat <- rep(NA, ntheta)
    theta_names <- c("gamma", "tau", "sigma", paste0("mu", 1:q), paste0("lambda", 1:nq))
    if(missing(Tz)) Tz <- Toeplitz(n = N)
    # profile likelihood on transformed scale
    ll.prof <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      tau <- itrans_tau(theta[2])
      sigma2 <- theta[3]^2
      acf1 <- fdyn_acf(alpha, tau, dT, N) + sigma2 * c(2, -1, rep(0, N-2))
      Tz$setAcf(acf1)
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      lmn.prof(suff)
    }
    # likelihood on transformed scale
    loglik <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      tau <- itrans_tau(theta[2])
      sigma2 <- theta[3]^2
      mu <- theta[3+1:q]
      Sigma <- itrans_Sigma(theta[3+q+1:nq]) # default: log(D)
      acf1 <- fdyn_acf(alpha, tau, dT, N) + sigma2 * c(2, -1, rep(0, N-2))
      Tz$setAcf(acf1)
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    }
    # calculate MLE
    fit <- optim(fn = ll.prof, par = c(-0.4, 0.1, 0.1),
                 control = list(fnscale = -1, ...))
    if(fit$convergence != 0) warning("optim did not converge.")
    theta_hat[1:3] <- fit$par # profiled parameters
    acf1 <- fdyn_acf(itrans_alpha(theta_hat[1]), itrans_tau(theta_hat[2]), dT, N) + 
      theta_hat[3]^2 * c(2, -1, rep(0, N-2))
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    theta_hat[3+1:q] <- suff$Beta
    theta_hat[3+q+1:nq] <- trans_Sigma(suff$S/suff$n)
    names(theta_hat) <- theta_names
    ans <- theta_hat # no-copy unless ans is modified
  }
  if(type == "dynamic") {
    ntheta <- 2+q+nq
    theta_hat <- rep(NA, ntheta)
    theta_names <- c("gamma", "tau", paste0("mu", 1:q), paste0("lambda", 1:nq))
    if(missing(Tz)) Tz <- Toeplitz(n = N)
    # profile likelihood on transformed scale
    ll.prof <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      tau <- itrans_tau(theta[2])
      acf1 <- fdyn_acf(alpha, tau, dT, N)
      Tz$setAcf(acf1)
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      lmn.prof(suff)
    }
    # likelihood on transformed scale
    loglik <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      tau <- itrans_tau(theta[2])
      mu <- theta[2+1:q]
      Sigma <- itrans_Sigma(theta[2+q+1:nq]) # default: log(D)
      acf1 <- fdyn_acf(alpha, tau, dT, N)
      Tz$setAcf(acf1)
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    }
    # calculate MLE
    fit <- optim(fn = ll.prof, par = c(0,0.1),
                 control = list(fnscale = -1, ...))
    if(fit$convergence != 0) warning("optim did not converge.")
    theta_hat[1:2] <- fit$par # profiled parameters
    acf1 <- fdyn_acf(itrans_alpha(theta_hat[1]), itrans_tau(theta_hat[2]), dT, N)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    theta_hat[2+1:q] <- suff$Beta
    theta_hat[2+q+1:nq] <- trans_Sigma(suff$S/suff$n)
    names(theta_hat) <- theta_names
    ans <- theta_hat # no-copy unless ans is modified
  }
  if(type == "localization") {
    ntheta <- 2+q+nq
    theta_hat <- rep(NA, ntheta)
    theta_names <- c("gamma", "sigma", paste0("mu", 1:q), paste0("lambda", 1:nq))
    if(missing(Tz)) Tz <- Toeplitz(n = N)
    # profile likelihood on transformed scale
    ll.prof <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      sigma2 <- theta[2]^2
      acf1 <- fbm_acf(alpha, dT, N) + sigma2 * c(2, -1, rep(0, N-2))
      Tz$setAcf(acf1)
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      lmn.prof(suff)  
    }
    # likelihood on transformed scale
    loglik <- function(theta) {
      alpha <- itrans_alpha(theta[1])
      sigma2 <- theta[2]^2
      mu <- theta[2+1:q]
      Sigma <- itrans_Sigma(theta[2+q+1:nq]) # default: log(D)
      acf1 <- fbm_acf(alpha, dT, N) + sigma2 * c(2, -1, rep(0, N-2))
      Tz$setAcf(acf1)
      suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
      lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    }
    # calculate MLE
    fit <- optim(fn = ll.prof, par = c(0, 0.1),
                 control = list(fnscale = -1, ...))
    if(fit$convergence != 0) warning("optim did not converge.")
    theta_hat[1:2] <- fit$par # profiled parameters
    acf1 <- fbm_acf(itrans_alpha(theta_hat[1]), dT, N) + 
      theta_hat[2] * c(2, -1, rep(0, N-2))
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    theta_hat[2+1:q] <- suff$Beta
    theta_hat[2+q+1:nq] <- trans_Sigma(suff$S/suff$n)
    names(theta_hat) <- theta_names
    ans <- theta_hat # no-copy unless ans is modified
  }
  if(var_calc) {
    # variance estimate
    V_hat <- hessian(loglik, x = theta_hat)
    V_hat <- solveV(-V_hat)
    colnames(V_hat) <- theta_names
    rownames(V_hat) <- theta_names
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}
