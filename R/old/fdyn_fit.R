#' Fit the fractional AR(1) model.
#'
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param tau Camera exposure time divided by `dT`
#' @param Tz Optional Toeplitz matrix for intermediate calculations.
#' @param var_calc If `TRUE`, also estimate variance matrix.
#' @param ... Additional `control` arguments to `stats::optim`.
#' @return Vector of coefficients and possibly variance matrix on the transformed scale (see Details).
#' @details The fractional AR(1) model has the form
#' \deqn{
#' \Delta X_n = (1-\rho) \Delta Z_n + \rho \Delta X_{n-1},
#' }
#' where \eqn{\Delta Z_n} are increments of a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by `trans(rho) = logit(1-rho/2)`, `trans(mu) = mu`, [trans_alpha()], and [trans_Sigma()].
#' @export
fdyn_fit <- function(dX, dT, tau, Tz, var_calc = TRUE) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  nq <- if(q == 1) 1 else 3
  ntheta <- 1+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", paste0("mu", 1:q), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on regular scale
  ll.prof <- function(theta) {
    alpha <- theta
    Tz$setAcf(fdyn_acf(alpha, tau, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.prof(suff)
  }
  # likelihood on transformed scale
  loglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    mu <- theta[1+1:q]
    Sigma <- itrans_Sigma(theta[1+q+1:nq]) # default: log(D)
    Tz$setAcf(fdyn_acf(alpha, tau, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  theta_hat[1] <- optimize(f = ll.prof,
                           interval = c(.01, 1.99), maximum = TRUE)$maximum
  Tz$setAcf(fdyn_acf(theta_hat[1], tau, dT, N)) # profiled parameters
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[1] <- trans_alpha(theta_hat[1]) # normalized scale
  theta_hat[1+1:q] <- suff$Beta
  theta_hat[1+q+1:nq] <- trans_Sigma(suff$S/suff$n)
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
