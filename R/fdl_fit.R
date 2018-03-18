#' Fit the fractional dynamic localization error model.
#'
#' @name fdl_fit
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param Tz Optional Toeplitz matrix for intermediate calculations.
#' @param var_calc If \code{TRUE}, also estimate variance matrix.
#' @param theta0 Length-3 vector of initial values for \code{(alpha, tau, sigma)}.  Default value is \code{(1, .1, .1)}.
#' @param ... Additional arguments to \code{\link[stats]{optim}}.
#' @return Vector of coefficients and possibly variance matrix on the transformed scale (see Details).
#' @details The fractional AR(1) model has the form
#' \deqn{
#' X_n = 1/\tau \int_0^\tau Z_(n+s) ds + \sigma e_{n},
#' }
#' where \eqn{Z_n} is a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed
#' scale defined by \code{trans(tau) = logit(tau)}, \code{trans(sigma) = log(sigma)}, \code{trans(mu) = mu}, \code{\link{trans_alpha}}, and
#' \code{\link{trans_Sigma}}.
#'
#' The functions \code{fdy_fit} and \code{flo_fit} fit pure dynamic error and pur localization error models, respectively.
#' @export
fdl_fit <- function(dX, dT, Tz, var_calc = TRUE, theta0, ...) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  nq <- if(q == 1) 1 else 3
  ntheta <- 3+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", "eta1", "eta2", paste0("mu", 1:q),
                   paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on transformed scale
  negll.prof <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    tau <- itrans_tau(theta[2])
    sigma2 <- exp(2*theta[3])
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.prof(suff)
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    tau <- itrans_tau(theta[2])
    sigma2 <- exp(2*theta[3])
    mu <- theta[3+1:q]
    Sigma <- itrans_Sigma(theta[3+q+1:nq]) # default: log(D)
    acf1 <- fdyn_acf(alpha, tau, dT, N) # + sigma2 * c(2, -1, rep(0, N-2))
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  if(missing(theta0)) theta0 <- c(1, .1, .1)
  tpar <- c(trans_alpha(theta0[1]), trans_tau(theta0[2]), log(theta0[3]))
  fit <- optim(fn = negll.prof, par = tpar, ...)
  if(fit$convergence != 0) stop("optim did not converge.")
  theta_hat[1:3] <- fit$par # profiled parameters
  acf1 <- fdyn_acf(alpha = itrans_alpha(theta_hat[1]),
                   tau = itrans_tau(theta_hat[2]), dT, N)
  acf1[1:2] <- acf1[1:2] + exp(2*theta_hat[3]) * c(2, -1)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[3+1:q] <- suff$Beta
  theta_hat[3+q+1:nq] <- trans_Sigma(suff$S/suff$n)
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

#' @rdname fdl_fit
#' @export
fdy_fit <- function(dX, dT, Tz, var_calc = TRUE, ...) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  nq <- if(q == 1) 1 else 3
  ntheta <- 2+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", "eta", paste0("mu", 1:q), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on transformed scale
  negll.prof <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    tau <- itrans_tau(theta[2])
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.prof(suff)
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    tau <- itrans_tau(theta[2])
    mu <- theta[2+1:q]
    Sigma <- itrans_Sigma(theta[2+q+1:nq]) # default: log(D)
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  fit <- optim(fn = negll.prof, par = c(0,0.1), ...)
  if(fit$convergence != 0) stop("optim did not converge.")
  theta_hat[1:2] <- fit$par # profiled parameters
  acf1 <- fdyn_acf(itrans_alpha(theta_hat[1]), itrans_tau(theta_hat[2]), dT, N)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[2+1:q] <- suff$Beta
  theta_hat[2+q+1:nq] <- trans_Sigma(suff$S/suff$n)
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

#' @rdname fdl_fit
#' @export
flo_fit <- function(dX, dT, Tz, var_calc = TRUE, ...) {
  # memory allocation
  N <- nrow(dX)
  q <- ncol(dX)
  nq <- if(q == 1) 1 else 3
  ntheta <- 2+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", "eta", paste0("mu", 1:q), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on transformed scale
  negll.prof <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    sigma2 <- exp(2*theta[2])
    acf1 <- fbm_acf(alpha, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.prof(suff)
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    sigma2 <- exp(2*theta[2])
    mu <- theta[2+1:q]
    Sigma <- itrans_Sigma(theta[2+q+1:nq]) # default: log(D)
    acf1 <- fbm_acf(alpha, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  fit <- optim(fn = ll.prof, par = c(0, 0.1),
               control = list(fnscale = -1, ...))
  if(fit$convergence != 0) stop("optim did not converge.")
  theta_hat[1:2] <- fit$par # profiled parameters
  acf1 <- fbm_acf(itrans_alpha(theta_hat[1]), dT, N)
  acf1[1:2] <- acf1[1:2] + exp(2*theta_hat[2]) * c(2, -1)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[2+1:q] <- suff$Beta
  theta_hat[2+q+1:nq] <- trans_Sigma(suff$S/suff$n)
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
