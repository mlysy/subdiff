#' Cochran's multivariate Q-test for heterogeneity.
#'
#' Under \code{H0}, \code{k} independent estimators are of a common \code{p}-diemensional parameter; their precision-weighted sum-of-squares follows asymptotically a \code{chi^2} distribution with \code{k-p} degrees of freedom.
#'
#' @param est A \code{k x p} matrix of independent estimates of a common parameter.
#' @param ve A \code{(p x p x k)}-dimensional array of corresponding variance estimates.
#' @return A list with elements:
#' \describe{
#'   \item{\code{Q}}{The value of the test statistic.}
#'   \item{\code{df}}{The degrees of freedom of the test.}
#'   \item{\code{pval}}{The p-value of the test.}
#'   \item{\code{est}}{The combined estimate of the parameter under \code{H0}.}
#'   \item{\code{ve}}{The variance estimate of the combined estimator.}
#' }
#' @export
cochranMQ <- function(est, ve) {
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
       est = est0, ve = ve0)
}
