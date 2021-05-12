#' Cochran's multivariate Q-test for heterogeneity.
#'
#' Suppose we have `k` independent estimators of `p`-dimensional vectors `theta_1, ..., theta_k`.  Under the null hypothesis `H0` that the `theta_i` are all the same, the Cochran Q-statistic is the precision-weighted sum-of-squares of the estimators.  Asymptotically it follows a chi-square distribution with `p*(k-1)` degrees of freedom.
#'
#' @param est A `k x p` matrix of independent estimates of a common parameter.
#' @param ve A `(p x p x k)`-dimensional array of corresponding variance estimates.
#' @return A list with elements:
#' \describe{
#'   \item{`Q`}{The value of the test statistic.}
#'   \item{`df`}{The degrees of freedom of the test.}
#'   \item{`pval`}{The p-value of the test.}
#'   \item{`coef`}{The combined estimate of the parameter under `H0`.}
#'   \item{`vcov`}{The variance estimate of the combined estimator.}
#' }
#' @export
cochran_mq <- function(est, ve) {
  k <- nrow(est)
  p <- ncol(est)
  # estimate of common parameter and its variance
  IV <- array(NA, dim = c(p,p+1,k))
  for(ii in 1:k) {
    # compute each ve^{-1} and ve^{-1} * est
    IV[,,ii] <- solveV(V = ve[,,ii], x = cbind(diag(p), est[ii,]))
  }
  ve0 <- rowSums(IV, dims = 2)
  ve0 <- solveV(ve0[,1:p], cbind(diag(p), ve0[,p+1]))
  est0 <- ve0[,p+1] # common parameter estimate
  ve0 <- ve0[,1:p] # common variance estimate
  # Q statistic
  Q <- 0
  for(ii in 1:k) {
    de <- est[ii,] - est0
    Q[1] <- Q + crossprod(de, IV[,1:p,ii] %*% de)
  }
  df0 <- p*(k-1) # is this correct???
  list(Q = Q,
       df = df0,
       pval = pchisq(Q, df = df0, lower.tail = FALSE),
       coef = est0, vcov = ve0)
}
