#' Fit the fractional AR(1) model.
#'
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @param theta0 Length-2 vector of initial values for `(alpha, rho)`.  Default value is `(1, 0)`.
#' @template args-dots_optim
#' @template ret-cov_vcov
#' @details The fractional AR(1) model has the form
#' \deqn{
#' \Delta X_n = (1-\rho) \Delta Z_n + \rho \Delta X_{n-1},
#' }
#' where \eqn{\Delta Z_n} are increments of a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by `trans(rho) = logit(1-rho/2)`, `trans(mu) = mu`, [trans_alpha()], and [trans_Sigma()].
#' @export
far_fit <- function(dX, dT, Tz, var_calc = TRUE, theta0, ...) {
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
    rho <- itrans_rho(theta[2])
    dY <- ar_resid(dX, rho, denom = FALSE)
    Tz$setAcf(fbm_acf(alpha, dT, N))
    suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
    nll <- -lmn.prof(suff)# + log(1-rho)*N*q
    nll
    # penalty on rho
    #nll + log1pe(theta[2]) + log1pe(-theta[2])
  }
  # likelihood on transformed scale
  negloglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    rho <- itrans_rho(theta[2])
    mu <- (1-rho) * theta[2+1:q]
    Sigma <- (1-rho)^2 * itrans_Sigma(theta[2+q+1:nq]) # default: log(D)
    dY <- ar_resid(dX, rho, denom = FALSE)
    Tz$setAcf(fbm_acf(alpha, dT, N))
    suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
    nll <- -lmn.loglik(Beta = t(mu),
                       Sigma = Sigma, suff = suff) # + log(1-rho)*N*q
    nll
    # penalty on rho
    ## nll + log1pe(theta[2]) + log1pe(-theta[2])
  }
  # calculate MLE
  if(missing(theta0)) theta0 <- c(1, 0)
  tpar <- c(trans_alpha(theta0[1]), trans_rho(theta0[2]))
  fit <- optim(fn = negll.prof, par = tpar, ...)
  if(fit$convergence != 0) {
    stop("optim did not converge (error code = ", fit$convergence, ").")
  }
  theta_hat[1:2] <- fit$par # profiled parameters
  Tz$setAcf(fbm_acf(itrans_alpha(theta_hat[1]), dT, N))
  dY <- ar_resid(dX, rho = itrans_rho(theta_hat[2]), denom = FALSE)
  suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
  rho_hat <- itrans_rho(theta_hat[2])
  theta_hat[2+1:q] <- suff$Beta/(1-rho_hat)
  theta_hat[2+q+1:nq] <- trans_Sigma(suff$S/suff$n/(1-rho_hat)^2)
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

# Nested function: calculates the residuals of an AR(1) model
# (1-rho) dZ_t = dX_t - rho dX_(t-1)
# seems this should be done on a different scale, with log(D) corrected later.
# denom = FALSE uses dZ_t = dX_t - rho dX_(t-1),
# with fBM+drift parameters mu and Sigma needing to be transformed to
# mu -> mu/(1-rho)   and Sigma -> Sigma/(1-rho)^2
ar_resid <- function(dX, rho, denom = TRUE) {
  N <- nrow(dX)
  dZ <- dX # no-copy (i.e., no memory allocation)
  if(denom) dZ <- dZ/(1-rho)
  dZ[-1,] <- dZ[-1,] - rho * dZ[-N,]
  ## dZ <- dX / (1-rho)
  ## dZ[-1, ] <- (dX[-1, ] - rho * dX[-N, ]) / (1-rho)
  ## dZ[-1, ] <- dX[-1, ]/(1-rho) - rho * dX[-N, ]/(1-rho)
  dZ
}


## far_fit <- function(dX, dT, Tz, var_calc = TRUE, theta0, ...) {
##   # memory allocation
##   N <- nrow(dX)
##   q <- ncol(dX)
##   nq <- if(q == 1) 1 else 3
##   ntheta <- 2+q+nq
##   theta_hat <- rep(NA, ntheta)
##   theta_names <- c("gamma", "eta", paste0("mu", 1:q), paste0("lambda", 1:nq))
##   if(missing(Tz)) Tz <- Toeplitz(n = N)
##   # profile likelihood on transformed scale
##   negll.prof <- function(theta) {
##     alpha <- itrans_alpha(theta[1])
##     rho <- itrans_rho(theta[2])
##     dY <- ar_resid(dX, rho)
##     Tz$setAcf(fbm_acf(alpha, dT, N))
##     suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
##     nll <- -lmn.prof(suff) + log(1-rho)*N*q
##     nll
##     # penalty on rho
##     #nll + log1pe(theta[2]) + log1pe(-theta[2])
##   }
##   # likelihood on transformed scale
##   negloglik <- function(theta) {
##     alpha <- itrans_alpha(theta[1])
##     rho <- itrans_rho(theta[2])
##     mu <- theta[2+1:q]
##     Sigma <- itrans_Sigma(theta[2+q+1:nq]) # default: log(D)
##     dY <- ar_resid(dX, rho)
##     Tz$setAcf(fbm_acf(alpha, dT, N))
##     suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
##     nll <- -lmn.loglik(Beta = t(mu),
##                        Sigma = Sigma, suff = suff) + log(1-rho)*N*q
##     nll
##     # penalty on rho
##     ## nll + log1pe(theta[2]) + log1pe(-theta[2])
##   }
##   # calculate MLE
##   if(missing(theta0)) theta0 <- c(1, 0)
##   tpar <- c(trans_alpha(theta0[1]), trans_rho(theta0[2]))
##   fit <- optim(fn = negll.prof, par = tpar, ...)
##   if(fit$convergence != 0) {
##     stop("optim did not converge (error code = ", fit$convergence, ").")
##   }
##   theta_hat[1:2] <- fit$par # profiled parameters
##   Tz$setAcf(fbm_acf(itrans_alpha(theta_hat[1]), dT, N))
##   dY <- ar_resid(dX, rho = itrans_rho(theta_hat[2]))
##   suff <- lmn.suff(Y = dY, X = dT, acf = Tz)
##   theta_hat[2+1:q] <- suff$Beta
##   theta_hat[2+q+1:nq] <- trans_Sigma(suff$S/suff$n)
##   names(theta_hat) <- theta_names
##   ans <- theta_hat # no-copy unless ans is modified
##   if(var_calc) {
##     # variance estimate
##     V_hat <- hessian(negloglik, x = theta_hat)
##     V_hat <- solveV(V_hat)
##     colnames(V_hat) <- theta_names
##     rownames(V_hat) <- theta_names
##     ans <- list(coef = theta_hat, vcov = V_hat)
##   }
##   ans
## }
