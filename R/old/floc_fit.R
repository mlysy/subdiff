#' Fit the fBM model with dynamic and static errors.
#'
#' @name floc_fit
#' @template args-dX
#' @template args-dT
#' @param tau Scalar between 0 and 1, indicating the percentage of time for which the camera shutter is open in the dynamic error model.  Estimated if missing. See Details.
#' @param sigma2 Magnitude of static error.  Estimated if missing. See Details.
#' @template args-Tz
#' @template args-var_calc
#' @param penalty logic, employ a small penalty on `(tau, sigma2)` when `TRUE`
#' @template ret-cov_vcov 
#' @details The fBM + dynamic and localization error (fdl) model has the form
#' \deqn{
#' X_n = 1/\tau \int_0^\tau Z_(n+s) ds + \sigma e_{n},
#' }
#' where \eqn{Z_n} is a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by `trans(tau) = logit(tau)`, `trans(sigma) = log(sigma)`, `trans(mu) = mu`, [trans_alpha()], and [trans_Sigma()].
#' When put `tau` = 0, this model becomes fractional localization error model.
#' When put `sigma2` = 0, this model becomes fractional dynamic error model.
#' When put `(tau, sigma2)` = 0, this model becomes fractional Brownian model.
#' 
#' @export
floc_fit <- function(dX, dT, Tz, var_calc = TRUE) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  
  theta_names <- c("alpha", "tau", "sigma2", paste0("mu", 1:qq),
                   paste0("lambda", 1:nq))
  # acf function on transformed scale
  acf_func <- function(theta) {
    floc_acf(itrans_alpha(theta[1]), itrans_tau(theta[2]), exp(2*theta[3]), dT, N)
  }
  
  # estimation
  ans <- Tz_fit(3, acf_func, dX, dT, Tz, var_calc, penalty = TRUE)
  
  # add names
  if(var_calc) {
    names(ans$coef) <- colnames(ans$vcov) <- 
      rownames(ans$vcov) <- theta_names
  } else {
    names(ans) <- theta_names
  }
  ans
}
