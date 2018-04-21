#' Fit the fBM model.
#'
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @template ret-cov_vcov
#' @details The fBM model is of the form
#' \deqn{
#' X_t = \mu t + \Sigma^{1/2} Z_t,
#' }
#' where \eqn{Z_t} consists of \code{q = 1,2} iid fBM processes with \eqn{MSD_Z(t) = t^\alpha}.  The MLE and variance estimate are calculated on the transformed scale defined by \code{trans(mu) = mu}, \code{\link{trans_alpha}}, and \code{\link{trans_Sigma}}.
#' @export
fbm_fit <- function(dX, dT, Tz, var_calc = TRUE) {
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
    Tz$setAcf(fbm_acf(alpha, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.prof(suff)
  }
  # likelihood on transformed scale
  loglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    mu <- theta[1+1:qq]
    Sigma <- itrans_Sigma(theta[1+qq+1:nq]) # default: log(D)
    Tz$setAcf(fbm_acf(alpha, dT, N))
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
