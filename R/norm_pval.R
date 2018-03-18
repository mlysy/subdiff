#' P-value calculation for various normality tests.
#'
#' @param x Vector of observations.
#' @return p-value of the normality test.
#' @name norm_pval
#' @details The following normality tests are currently implemented:
#' \describe{
#'   \item{\code{ad_pval}}{The Anderson-Darling test, which assumes iid observations.  This is merely a wrapper to \code{\link[nortest]{ad.test}}.}
#'   \item{\code{sw_pval}}{A modified Shapiro-Wilk test, which also assumes iid observations.  For sample sizes less than 5000, this function is merely a wrapper to \code{\link[stats]{shapiro.test}}.  For sample sizes greater than 5000, the algorithm due to Royston (1992) used in \code{\link[stats]{shapiro.test}} breaks down.  The modification proposed here is to (i) randomly divide the sample into \code{K} roughly equal parts each less than 5000, (ii) calculate the Shapiro-Wilk p-value for each, and (iii) note that the minimum of \code{K} independent p-values has a \code{Beta(1, K)} distribution.  It is this second-stage p-value whic is returned by \code{sw_pval}.}
#'   \item{\code{bc_pval}}{The so-called "Berkowitz' Correlation" test.  That is, An autoregressive model of the form
#' \preformatted{
#' x[n] = rho * x[n-1] + z[n],    z[n] ~iid N(0,1)
#' }
#' is fit to the data, and a likelihood ratio test against \code{H0: rho = 0} is performed.}
#' }
#' @note The tests above **only** check for normality, not that \code{x} is iid \code{N(0,1)}.  Such modifications will need to be made before the package is released.
#' @references Royston, P. "Approximating the Shapiro-Wilk W-test for non-normality." \emph{Statistics and Computing} 2:3 (1992): 117-119. \url{https://doi.org/10.1007/BF01891203}.
#'
#' Berkowitz, J. "Testing density forecasts, with applications to risk management." \emph{Journal of Business and Economic Statistics} 19:4 (2001): 465-474. \url{https://doi.org/10.1198/07350010152596718}.

# anderson-darling test
#' @rdname norm_pval
#' @export
ad_pval <- function(x) ad.test(x)$p.value

# shapiro-wilk test
#' @rdname norm_pval
#' @export
sw_pval <- function(x) {
  nmax <- 5000 # maximum test size.  divide into subtests when n > 5000
  n <- length(x)
  k <- ceiling(n/nmax)
  ind <- sample(rep(1:k, each = ceiling(n/k), len = n))
  pv <- tapply(x, ind, function(z) shapiro.test(z)$p.value)
  pbeta(min(pv), shape1 = 1, shape2 = length(pv))
}


# berkowitz-correlation test
#' @rdname norm_pval
#' @export
bc_pval <- function(x) {
  n <- length(x)-1
  x1 <- x[1:n]
  x2 <- x[1+1:n]
  # regression coefficients
  b1 <- cov(x2, x1)/var(x1)
  b00 <- mean(x2)
  b0 <- b00 - b1 * mean(x1)
  # residuals
  z2 <- (x2 - b0 - b1*x1)
  z1 <- (x2 - b00)
  # likelihood ratio statistic
  lambda <- n * (log(mean(z1^2)) - log(mean(z2^2)))
  pchisq(lambda, df = 1) # p-value
}
