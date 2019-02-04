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
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  theta_names <- c("alpha", paste0("mu", 1:qq), paste0("lambda", 1:nq))
  
  # acf function on transformed scale
  acf_func <- function(theta) {
    fbm_acf(itrans_alpha(theta), dT, N)
  }
  
  # estimation
  ans <- Tz_fit(1, acf_func, dX, dT, Tz, var_calc)

  # add names
  if(var_calc) {
    names(ans$coef) <- colnames(ans$vcov) <- 
      rownames(ans$vcov) <- theta_names
  } else {
    names(ans) <- theta_names
  }
  ans
}
