#' Fit the fractional MA(1) model.
#'
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param nlag Maximum number of lags for the fMA model.
#' @param Tz Optional Toeplitz matrix for intermediate calculations.
#' @param var_calc If \code{TRUE}, also estimate variance matrix.
#' @param ... Additional arguments to \code{stats::optim}.
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
  negll.prof <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    rho <- itrans_rho(theta[1+1:nlag])
    Tz$setAcf(fma_acf(alpha, rho, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.prof(suff)
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    rho <- itrans_rho(theta[1+1:nlag])
    mu <- theta[1+nlag+1:q]
    Sigma <- itrans_Sigma(theta[1+nlag+q+1:nq]) # default: log(D)
    Tz$setAcf(fma_acf(alpha, rho, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  fit <- optim(fn = negll.prof, par = rep(0,nlag+1), ...)
  if(fit$convergence != 0) stop("optim did not converge.")
  theta_hat[1+0:nlag] <- fit$par # profiled parameters
  Tz$setAcf(fma_acf(itrans_alpha(theta_hat[1]), itrans_rho(theta_hat[1+1:nlag]), dT, N))
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[1+nlag+1:q] <- suff$Beta
  theta_hat[1+nlag+q+1:nq] <- trans_Sigma(suff$S/suff$n)
  names(theta_hat) <- theta_names
  ans <- theta_hat # no-copy unless ans is modified
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


#' @rdname fma_fit
#' @export
fma_trunc <- function(dX, dT, rho) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  theta_hat <- rep(NA, 1)
  theta_names <- c("gamma")
  Tz <- Toeplitz(n = N)
  # profile likelihood on regular scale
  ll.prof <- function(theta) {
    alpha <- itrans_alpha(theta)
    Tz$setAcf(fma_acf(alpha, rho, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.prof(suff)
  }
  # calculate MLE
  theta_hat <- optimize(f = ll.prof,
                        interval = c(-5.3, -0.02), maximum = TRUE)$maximum
  names(theta_hat) <- theta_names
  # variance estimate
  V_hat <- hessian(ll.prof, x = theta_hat)
  V_hat <- subdiff:::solveV(-V_hat)
  colnames(V_hat) <- theta_names
  rownames(V_hat) <- theta_names
  list(coef = theta_hat, vcov = V_hat)
}
