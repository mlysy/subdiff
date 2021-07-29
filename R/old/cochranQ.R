#' Cochran's Q-test for heterogeneity.
#'
#' Under `H0`, `k` independent estimators are of a common parameter; their precision-weighted sum-of-squares follows asymptotically a chi^2 distribution with `k-1` degrees of freedom.
#'
#' @param est A vector of independent estimates of a common parameter.
#' @param se A corresponding vector of standard errors.
#' @return The value of the test statistic, its p-value, and the weighted mean.
#' @export
cochranQ <- function(est, se) {
  k <- length(est)
  wgt <- 1/se^2
  ebar <- sum(wgt * est)/sum(wgt)
  Q <- sum(wgt*(est - ebar)^2)
  c(Q = Q, pval = pchisq(Q, df = k-1, lower.tail = FALSE), ebar = ebar)
}
