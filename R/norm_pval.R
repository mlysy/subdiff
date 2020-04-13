#' P-value calculation for various normality tests.
#'
#' @param x Vector of observations.
#' @param stdn Logical, whether or not to test against standard normal distribution, or take mean and variance to be unknown.
#' @return p-value of the normality test.
#' @name norm_pval
#' @details The following normality tests are currently implemented:
#' \describe{
#'   \item{`ad_pval`}{The Anderson-Darling test, which assumes iid observations.  This is merely a wrapper to `nortest::ad.test`.}
#'   \item{`sw_pval`}{A modified Shapiro-Wilk test, which also assumes iid observations.  For sample sizes less than 5000, this function is merely a wrapper to [stats::shapiro.test()].  For sample sizes greater than 5000, the algorithm due to Royston (1992) used in `stats::shapiro.test` breaks down.  The modification proposed here is to (i) randomly divide the sample into `K` roughly equal parts each less than 5000, (ii) calculate the Shapiro-Wilk p-value for each, and (iii) note that the minimum of `K` independent p-values has a `Beta(1, K)` distribution.  It is this second-stage p-value whic is returned by `sw_pval`.}
#'   \item{`bc_pval`}{The so-called "Berkowitz' Correlation" test.  That is, An autoregressive model of the form
#' \preformatted{
#' x[n] = rho * x[n-1] + z[n],    z[n] ~iid N(0,1)
#' }
#' is fit to the data, and a likelihood ratio test against `H0: rho = 0` is performed.}
#' }
#' @note The tests above **only** check for normality, not that `x` is iid `N(0,1)`.  Such modifications will need to be made before the package is released.
#' @references Royston, P. "Approximating the Shapiro-Wilk W-test for non-normality." *Statistics and Computing* 2:3 (1992): 117-119. <https://doi.org/10.1007/BF01891203>.
#'
#' Berkowitz, J. "Testing density forecasts, with applications to risk management." *Journal of Business and Economic Statistics* 19:4 (2001): 465-474. <https://doi.org/10.1198/07350010152596718>.

# anderson-darling test
#' @rdname norm_pval
#' @export
ad_pval <- function(x, stdn = TRUE) {
  x <- sort(x)
  n <- length(x)
  if (n < 8) stop("sample size must be greater than 7")
  if(stdn) {
    mu <- 0
    sig <- 1
  } else {
    mu <- mean(x)
    sig <- sd(x)
  }
  logp1 <- pnorm((x - mu)/sig, log.p = TRUE)
  logp2 <- pnorm(-(x - mu)/sig, log.p = TRUE)
  h <- (2 * seq(1:n) - 1) * (logp1 + rev(logp2))
  A <- -n - mean(h)
  AA <- (1 + 0.75/n + 2.25/n^2) * A
  if(AA < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
  } else if(AA < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
  } else if(AA < 0.6) {
    pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
  } else if(AA < 10) {
    pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
  } else pval <- 3.7e-24
  pval
}

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
