#' Fit the fractional dynamic localization error model.
#'
#' @name fdl_fit
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param tau Plug-in coefficient, requires estimation if missing,
#' @param sigma2 Plug-in coefficient, requires estimation if missing,
#' @param Tz Optional Toeplitz matrix for intermediate calculations.
#' @param var_calc If \code{TRUE}, also estimate variance matrix.
#' @param penalty logic, employ a small penalty on \code{(tau, sigma2)} when \code{TRUE}
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
#' When put \code{tau} = 0, this model becomes fractional localization error model.
#' When put \code{sigma2} = 0, this model becomes fractional dynamic error model.
#' When put \code{(tau, sigma2)} = 0, this model becomes fractional Brownian model.
#' 
#' @export
fdl_fit <- function(dX, dT, tau, sigma2, Tz, var_calc = TRUE, penalty = TRUE, theta0, ...) {
  if(missing(tau) & missing(sigma2)) {
    ans <- fdylo_fit(dX, dT, Tz, var_calc, penalty, theta0, ...)
  } else if(missing(tau)) {
    ans <- fdy_fit(dX, dT, sigma2, Tz, var_calc, penalty, ...)
  } else if(missing(sigma2)) {
    ans <- flo_fit(dX, dT, tau, Tz, var_calc, penalty, ...)
  } else {
    ans <- f_fit(dX, dT, tau, sigma2, Tz, var_calc)
  }
  ans
}

fdylo_fit <- function(dX, dT, Tz, var_calc, penalty, theta0, ...) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  ntheta <- 3+qq+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", "tau", "sigma2", paste0("mu", 1:qq),
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
    nlp <- -lmn.prof(suff)
    if(penalty) {
      # penalty on tau, sigma
      nlp <- nlp + log1pe(theta[2]) + log1pe(-theta[2]) - theta[3]
    }
    nlp
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    tau <- itrans_tau(theta[2])
    sigma2 <- exp(2*theta[3])
    mu <- theta[3+1:qq]
    Sigma <- itrans_Sigma(theta[3+qq+1:nq]) # default: log(D)
    acf1 <- fdyn_acf(alpha, tau, dT, N) # + sigma2 * c(2, -1, rep(0, N-2))
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    nlp <- -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    if(penalty) {
      # penalty on tau, sigma
      nlp <- nlp + log1pe(theta[2]) + log1pe(-theta[2]) - theta[3]
    }
    nlp
  }
  # calculate MLE
  if(missing(theta0)) theta0 <- c(1, .1, .1)
  tpar <- c(trans_alpha(theta0[1]), trans_tau(theta0[2]), log(theta0[3]))
  fit <- optim(fn = negll.prof, par = tpar, ...)
  if(fit$convergence != 0) {
    stop("optim did not converge (error code = ", fit$convergence, ").")
  }
  theta_hat[1:3] <- fit$par # profiled parameters
  acf1 <- fdyn_acf(alpha = itrans_alpha(theta_hat[1]),
                   tau = itrans_tau(theta_hat[2]), dT, N)
  acf1[1:2] <- acf1[1:2] + exp(2*theta_hat[3]) * c(2, -1)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[3+1:qq] <- suff$Beta
  theta_hat[3+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
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

fdy_fit <- function(dX, dT, sigma2, Tz, var_calc, penalty, theta0, ...) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  ntheta <- 2+qq+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", "tau", paste0("mu", 1:qq), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on transformed scale
  negll.prof <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    tau <- itrans_tau(theta[2])
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    nlp <- -lmn.prof(suff)
    if(penalty) {
      # penalty on tau
      nlp <- nlp + log1pe(theta[2]) + log1pe(-theta[2])
    }
    nlp
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    tau <- itrans_tau(theta[2])
    mu <- theta[2+1:qq]
    Sigma <- itrans_Sigma(theta[2+qq+1:nq]) # default: log(D)
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    nlp <- -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    if(penalty) {
      # penalty on tau
      nlp <- nlp + log1pe(theta[2]) + log1pe(-theta[2])
    }
    nlp
  }
  # calculate MLE
  if(missing(theta0)) theta0 <- c(1, .1)
  tpar <- c(trans_alpha(theta0[1]), trans_tau(theta0[2]))
  fit <- optim(fn = negll.prof, par = tpar, ...)
  if(fit$convergence != 0) stop("optim did not converge.")
  theta_hat[1:2] <- fit$par # profiled parameters
  acf1 <- fdyn_acf(itrans_alpha(theta_hat[1]), itrans_tau(theta_hat[2]), dT, N)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[2+1:qq] <- suff$Beta
  theta_hat[2+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
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

flo_fit <- function(dX, dT, tau, Tz, var_calc, penalty, theta0, ...) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  ntheta <- 2+qq+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", "sigma2", paste0("mu", 1:qq), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on transformed scale
  negll.prof <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    sigma2 <- exp(2*theta[2])
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    nlp <- -lmn.prof(suff)
    if(penalty) {
      # penalty on sigma
      nlp <- nlp - theta[3]
    }
    nlp
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    sigma2 <- exp(2*theta[2])
    mu <- theta[2+1:qq]
    Sigma <- itrans_Sigma(theta[2+qq+1:nq]) # default: log(D)
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    nlp <- -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
    if(penalty) {
      # penalty on sigma
      nlp <- nlp - theta[3]
    }
    nlp
  }
  # calculate MLE
  if(missing(theta0)) theta0 <- c(1, .1)
  tpar <- c(trans_alpha(theta0[1]), log(theta0[3]))
  fit <- optim(fn = negll.prof, par = c(0,.1),
               control = list(fnscale = -1, ...))
  if(fit$convergence != 0) stop("optim did not converge.")
  theta_hat[1:2] <- fit$par # profiled parameters
  acf1 <- fbm_acf(itrans_alpha(theta_hat[1]), dT, N)
  acf1[1:2] <- acf1[1:2] + exp(2*theta_hat[2]) * c(2, -1)
  Tz$setAcf(acf1)
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[2+1:qq] <- suff$Beta
  theta_hat[2+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
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

f_fit <- function(dX, dT, tau, sigma2, Tz, var_calc) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  ntheta <- 1+qq+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", paste0("mu", 1:qq), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on regular scale
  ll.prof <- function(theta) {
    alpha <- theta
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.prof(suff)
  }
  # likelihood on transformed scale
  loglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    acf1 <- fdyn_acf(alpha, tau, dT, N)
    acf1[1:2] <- acf1[1:2] + sigma2 * c(2, -1)
    mu <- theta[1+1:qq]
    Sigma <- itrans_Sigma(theta[1+qq+1:nq]) # default: log(D)
    Tz$setAcf(acf1)
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  theta_hat[1] <- optimize(f = ll.prof,
                           interval = c(.01, 1.99), maximum = TRUE)$maximum
  Tz$setAcf(fbm_acf(theta_hat[1], dT, N)) # profiled parameters
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[1] <- trans_alpha(theta_hat[1]) # normalized scale
  theta_hat[1+1:qq] <- suff$Beta
  theta_hat[1+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
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