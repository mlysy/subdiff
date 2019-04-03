#' Fit the fractional Autoregressive Moving-average model.
#'
#' @template args-dX
#' @template args-dT
#' @param nlag Number of lags (see Details).
#' @param alpha Subdiffusion parameter, optional.
#' @template args-Tz
#' @template args-var_calc
#' @template args-dots_optim
#' @template ret-cov_vcov
#' @details The fractional Moving-average model has the form
#' \deqn{
#' \Delta X_n = (1-\sum_{i=1}^p \rho_i) \Delta Z_n + \sum_{i=1}^p \rho_{i} \Delta Z_{n-i},
#' }
#' where \eqn{\Delta Z_n} are increments of a 1D or 2D fBM process. The MLE and variance estimate are calculated on the transformed scale defined by \code{trans(rho) = logit(1-rho/2)}, \code{trans(mu) = mu}, \code{\link{trans_alpha}}, and \code{\link{trans_Sigma}}.
#' @export
farma_fit <- function(dX, dT, order, Tz, var_calc = TRUE) {
  # memory allocation
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- if(qq == 1) 1 else 3
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  
  if(order[1]) {
    # order = c(p,q): arma
    theta_names <- c("alpha", paste0("theta", 1:order[1]), 
                     paste0("rho", 1:order[2]), paste0("mu", 1:qq),
                     paste0("lambda", 1:nq))
    # acf function on transformed scale
    acf_func <- function(theta) {
      farma_acf(itrans_alpha(theta[1]), itrans_rho(theta[2]), itrans_rho(theta[2+1:order[2]]), dT, N)
    }
  } else {
    # order = c(0,q): ma
    theta_names <- c("alpha", paste0("rho", 1:order[2]), paste0("mu", 1:qq),
                     paste0("lambda", 1:nq))
    # acf function on transformed scale
    acf_func <- function(theta) {
      farma_acf(itrans_alpha(theta[1]), 0, itrans_rho(theta[1+1:order[2]]), dT, N)
    }
  }
  
  # estimation
  ans <- Tz_fit(1+sum(order), acf_func, dX, dT, Tz, var_calc)
  
  # add names
  if(var_calc) {
    names(ans$coef) <- colnames(ans$vcov) <- 
      rownames(ans$vcov) <- theta_names
  } else {
    names(ans) <- theta_names
  }
  ans
}

