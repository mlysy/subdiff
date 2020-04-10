#' Fit the location-scale model for Gaussian CSI process.
#'
#' Estimate the coefficients and their covariance for parameters in location-scale model where target data follows Gaussian CSI process (See \strong{Details}).
#'
#' @param model An list of class \code{csi_class} (see \code{\link{fbm_model}}, \code{\link{floc_model}}, \code{\link{farma_model}}).
#' @template args-dX
#' @template args-dT
#' @template args-Tz
#' @template args-var_calc
#' @template ret-cov_vcov
#'
#' @details The location-scale model is of following form
#' \deqn{
#' X_n = \mu n \Delta t + \Sigma^{1/2} Z_n
#' }{
#' X[n] = \mu n \Delta t + \Sigma^{1/2} Z[n]
#' }
#' where \eqn{\mu} is the drift parameter, \eqn{\Sigma} is the between-trajectory covariance and \eqn{Z[n]} is a Gaussian continuous-stationary-increment (CSI) process.
#'
#' In the location-scale model, \eqn{\mu, \Sigma} are nuisance parameters and can be profiled out using function \code{LMN::lmn_prof}.
#'
#' @references Lysy, M., Pillai, N.S., Hill, D.B., Forest, M.G., Mellnik, J.W.R., Vasquez, P.A., and McKinley, S.A. "Model comparison and assessment for single particle tracking in biological fluids." \emph{Journal of the American Statistical Association} 111.516 (2016): 1413-1426. \url{https://doi.org/10.1080/01621459.2016.1158716}.
#'
#' @example examples/fit_setup.R
#' @example examples/csi_fit.R
#'
#' @export
csi_fit <- function(model, dX, dT, Tz, var_calc) {
  # problem dimensions
  N <- nrow(dX)
  qq <- ncol(dX)
  nq <- get_nq(qq)
  if(missing(Tz)) Tz <- Toeplitz$new(N = N)

  # parameter setup
  ntheta <- length(model$theta_names)
  theta_hat <- rep(NA, ntheta+qq+nq)
  tnames <- c(model$theta_names,
                        paste0("mu", 1:qq),
                        paste0("lambda", 1:nq))
  names(theta_hat) <- tnames

  # profile likelihood on transformed scale
  negll_prof <- function(gamma) {
    # convert theta into original scale
    theta <- model$theta_itrans(gamma)
    # ACF
    acf1 <- model$acf(theta, dT, N)
    Tz$set_acf(acf1)
    # profile likelihood
    suff <- lmn_suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
    nlp <- -lmn_prof(suff)
    # penalty
    nlp <- nlp + model$penalty(gamma)
    nlp
  }

  # calculate MLE on transformed scale.
  if(ntheta == 1) {
    theta_hat[1:ntheta] <- optimize(f = negll_prof, interval = c(-5, 5))$minimum
  } else {
    fit <- optim(fn = negll_prof, par = rep(0,ntheta))
    if(fit$convergence != 0) warning("optim did not converge.")
    theta_hat[1:ntheta] <- fit$par
  }

  # ACF
  acf1 <- model$acf(model$theta_itrans(theta_hat[1:ntheta]), dT, N)
  Tz$set_acf(acf1)

  # profile likelihood
  suff <- lmn_suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
  theta_hat[ntheta+1:qq] <- suff$Bhat
  theta_hat[ntheta+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
  ans <- theta_hat

  if(var_calc) {
    # likelihood on transformed scale
    negloglik <- function(gamma) {
      theta <- model$theta_itrans(gamma[1:ntheta])
      mu <- gamma[ntheta+1:qq]
      Sigma <- itrans_Sigma(gamma[ntheta+qq+1:nq])
      acf1 <- model$acf(theta, dT, N)
      Tz$set_acf(acf1)
      suff <- lmn_suff(Y = dX, X = dT, V = Tz, Vtype = "acf")
      nlp <- -lmn_loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
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

