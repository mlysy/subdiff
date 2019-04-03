#' Fitting location-scale model for Gaussian CSI process.
#' 
#' @name csi_fit
#' @param model An S3 object of type `csi_class` (see \code{\link{fbm_model}}, \code{\link{floc_model}}, \code{\link{farma_model}})
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @template ret-cov_vcov
#' 
#' @details The location-scale model for CSI is of the form
#' \deqn{
#' \Delta X_n = \mu \Delta t + \Sigma^{1/2} \Delta Z_n,
#' }
#' where \eqn{\Delta Z_n = Z_n - Z_{n-1}} is the increment process of CSI process.
#' 
#' @example example/parameters.R
#' @examples 
#' model <- fbm_model() # generating fBM model object
#' fits <- csi_fit(model, dX, dT, Tz, var_calc = TRUE)
#' 
#' @export
csi_fit <- function(model, dX, dT, Tz, var_calc) {
  # problem dimensions
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- get_nq(qq)
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  
  # parameter setup
  ntheta <- length(model$theta_names)
  theta_hat <- rep(NA, ntheta+qq+nq)
  tnames <- c(model$theta_names, 
                        paste0("mu", 1:qq), 
                        paste0("lambda", 1:nq))
  names(theta_hat) <- tnames
  
  # profile likelihood on transformed scale
  negll.prof <- function(gamma) {
    # convert theta into original scale
    theta <- model$theta_itrans(gamma)
    # ACF
    acf1 <- model$acf(theta, dT, N)
    Tz$setAcf(acf1)
    # profile likelihood
    suff <- lmn.suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
    nlp <- -lmn.prof(suff)
    # penalty
    nlp <- nlp + model$penalty(gamma)
    nlp
  }

  # calculate MLE on transformed scale.
  if(ntheta == 1) {
    theta_hat[1:ntheta] <- optimize(f = negll.prof, interval = c(-5, 5))$minimum
  } else {
    fit <- optim(fn = negll.prof, par = rep(0,ntheta))
    if(fit$convergence != 0) warning("optim did not converge.")
    theta_hat[1:ntheta] <- fit$par
  }
  
  # ACF
  acf1 <- model$acf(model$theta_itrans(theta_hat[1:ntheta]), dT, N)
  Tz$setAcf(acf1)
  
  # profile likelihood
  suff <- lmn.suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
  theta_hat[ntheta+1:qq] <- suff$Beta
  theta_hat[ntheta+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
  ans <- theta_hat
  
  if(var_calc) {
    # likelihood on transformed scale
    negloglik <- function(gamma) {
      theta <- model$theta_itrans(gamma[1:ntheta])
      mu <- gamma[ntheta+1:qq]
      Sigma <- itrans_Sigma(gamma[ntheta+qq+1:nq])
      acf1 <- model$acf(theta, dT, N)
      Tz$setAcf(acf1)
      suff <- lmn.suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
      nlp <- -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
      nlp
    }
    # variance estimate
    V_hat <- hessian(negloglik, x = theta_hat)
    V_hat <- solveV(V_hat)
    colnames(V_hat) <- rownames(V_hat) <- tnames
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}

