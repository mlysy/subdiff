
#--- notes -----------------------------------------------------------------

# need:

# 1. boot_conf
# 2. msd_conf

#' Bootstrap standard errors and confidence intervals.
#'
#' @param mu_boot Matrix of bootstrap samples; each row is a bootstrap estimate of a parameter vector \code{mu}.
#' @param mu_hat Vector of point estimates for each element of \code{mu}.  Defaults to \code{colMeans(mu_boot)}.
#' @param conf Percent size of confidence interval (scalar between 0 and 1).
#' @return Matrix with \code{length(mu)} rows and columns named:
#' \describe{
#'   \item{\code{coef}}{Point estimates of \code{mu}.}
#'   \item{\code{se}}{Standard errors.}
#'   \item{\code{L,U}}{Endpoints of confidence intervals.}
#' }
#' @export
boot_conf <- function(mu_boot, mu_hat = NULL, conf = .95) {
  conf <- (1-conf)/2
  if(is.null(mu_hat)) mu_hat <- colMeans(mu_boot)
  se <- apply(mu_boot, 2, sd)
  ci <- apply(2*mu_hat - t(mu_boot), 1, quantile, probs = c(conf, 1-conf))
  cbind(coef = mu_hat, se = se, L = ci[1,], U = ci[2,])
}
