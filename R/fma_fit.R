#' Fit the fractional MA(1) model.
#'
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param Tz Optional Toeplitz matrix for intermediate calculations.
#' @param var_calc If \code{TRUE}, also estimate variance matrix.
#' @param ... Additional \code{control} arguments to \code{stats::optim}.
#' @return Vector of coefficients and possibly variance matrix on the transformed scale (see Details).
#' @details The fractional MA(1) model has the form
#' \deqn{
#' \Delta X_n = (1-\rho) \Delta Z_n + \rho \Delta Z_{n-1},
#' }
#' where \eqn{\Delta Z_n} are increments of a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by \code{trans(rho) = logit(1-rho/2)}, \code{trans(mu) = mu}, \code{\link{trans_alpha}}, and \code{\link{trans_Sigma}}.
#' @export
fma_fit <- function(dX, dT, nlag, Tz, var_calc = TRUE, ...) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  nq <- if(q == 1) 1 else 3
  if(missing(nlag)) nlag <- 1
  ntheta <- 1+nlag+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", paste0("eta", 1:nlag), paste0("mu", 1:q), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on transformed scale
  ll.prof <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    rho <- itrans_rho(theta[1+1:nlag])
    dY <- ma_resid(dX, rho)
    Tz$setAcf(fbm_acf(alpha, dT, N))
    suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
    lmn.prof(suff) - log(1-sum(rho))*N*q
  }
  # likelihood on transformed scale
  loglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    rho <- itrans_rho(theta[1+1:nlag])
    mu <- theta[1+nlag+1:q]
    Sigma <- itrans_Sigma(theta[1+nlag+q+1:nq]) # default: log(D)
    dY <- ma_resid(dX, rho)
    Tz$setAcf(fbm_acf(alpha, dT, N))
    suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
    lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff) - log(1-sum(rho))*N*q
  }
  # calculate MLE
  fit <- optim(fn = ll.prof, par = rep(0,nlag+1),
               control = list(fnscale = -1, ...))
  if(fit$convergence != 0) warning("optim did not converge.")
  theta_hat[1+0:nlag] <- fit$par # profiled parameters
  Tz$setAcf(fbm_acf(itrans_alpha(theta_hat[1]), dT, N))
  dY <- ma_resid(dX, rho = itrans_rho(theta_hat[1+1:nlag]))
  suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
  theta_hat[1+nlag+1:q] <- suff$Beta
  theta_hat[1+nlag+q+1:nq] <- trans_Sigma(suff$S/suff$n)
  names(theta_hat) <- theta_names
  ans <- theta_hat # no-copy unless ans is modified
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

#' Nested function
ma_resid <- function(dX, rho) {
  if(length(rho) == 1) {
    ans <- ma1_resid(dX, rho)
  } else if(length(rho) == 2) {
    ans <- ma2_resid(dX, rho[1], rho[2])
  } else if(length(rho) == 3) {
    ans <- ma3_resid(dX, rho[1], rho[2], rho[3])
  } else{
    stop("Lag of fMA model is at most 3")
  }
  ans
}