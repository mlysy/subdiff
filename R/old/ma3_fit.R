#' Fit the MA(3) model.
#'
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param var_calc If `TRUE`, also estimate variance matrix.
#' @param ... Additional `control` arguments to `stats::optim`.
#' @return Vector of coefficients and possibly variance matrix on the transformed scale (see Details).
#' @details The fractional MA(3) model has the form
#' \deqn{
#' \Delta X_n =  \mu \Delta t + \Sigma^{1/2} Z_n,
#' }
#' \deqn{
#' Z_n = \epsilon_n + \phi_1 \epsilon_{n-1} + \phi_2 \epsilon_{n-2} + \phi_3 \epsilon_{n-3},
#' }
#' where \eqn{\epsilon_n \sim N(0,1)}. The MLE and variance estimate are calculated on the transformed scale defined by `trans(phi) = phi`, `trans(mu) = mu`, and [trans_Sigma()].  Fitting is done by assuming that `eps[0] = eps[-1] = eps[-2] = 0`, leading to a simple formula for the model residuals via the `stats::filter` function.  Preliminary parameter estimates obtained from `stats::arima` are passed on to `optim`.
#' @export
ma3_fit <- function(dX, dT, var_calc = TRUE, ...) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  nq <- if(q == 1) 1 else 3
  ntheta <- 3+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c(paste0("phi", 1:3),
                   paste0("mu", 1:q), paste0("lambda", 1:nq))
  # profile likelihood on transformed scale
  ll.prof <- function(theta) {
    phi <- theta[1:3]
    suff <- ma3_suff(dX, dT, phi)
    lmn.prof(suff)
  }
  # likelihood on transformed scale
  loglik <- function(theta) {
    phi <- theta[1:3]
    mu <- theta[3+1:q]
    Sigma <- itrans_Sigma(theta[3+q+1:nq]) # default: log(D)
    suff <- ma3_suff(dX, dT, phi)
    lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  fit <- arima(dX[,1], order = c(0,0,3), method = "CSS") # preliminary fit
  fit <- optim(fn = ll.prof, par = coef(fit)[1:3],
               control = list(fnscale = -1, ...))
  if(fit$convergence != 0) stop("optim did not converge.")
  theta_hat[1:3] <- fit$par # profiled parameters
  suff <- ma3_suff(dX, dT, theta_hat[1:3])
  theta_hat[3+1:q] <- suff$Beta
  theta_hat[3+q+1:nq] <- trans_Sigma(suff$S/suff$n)
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

# lmn-type sufficient statistics for ma3 model
# should modify lmn to accept IP/ldV
ma3_suff <- function(Y, dT, phi) {
  n <- nrow(Y)
  q <- ncol(Y)
  npred <- 0
  X <- matrix(dT, n, 1)
  noBeta <- FALSE
  p <- (!noBeta) * ncol(X)
  # inner product calculations
  Z <- matrix(0, n, q+p+npred)
  Z[,1:q] <- Y
  if(!noBeta) {
    Z[,q+(1:p)] <- X[1:n,]
  }
  # inner product calculations
  IP <- crossprod(ma3_solve(Z, phi))
  ldV <- 0
  # convert inner products to sufficient statistics
  S <- IP[1:q,1:q,drop=FALSE]
  if(noBeta) {
    Beta.hat <- NULL
    T <- NULL
  } else {
    T <- IP[q+(1:p),q+(1:p),drop=FALSE]
    Beta.hat <- solve(T, IP[q+(1:p),1:q,drop=FALSE])
    S <- S - IP[1:q,q+(1:p)] %*% Beta.hat
  }
  # output
  list(Beta.hat = Beta.hat, S = S, T = T, ldV = ldV, n = n)
}

# invert the ma(3) model
ma3_solve <- function(X, theta) {
  # problem dimensions
  X <- as.matrix(X)
  dd <- ncol(X)
  N <- nrow(X)
  eps <- filter(X, filter = -theta, method = "recursive")
  tsp(eps) <- NULL
  eps
}
