#' Fit the fSD model.
#'
#' @template args-Xt
#' @template args-dt
#' @template args-drift_preset
#' @template args-vcov
#' @template args-ad_only
#'
#' @template ret-fit
#'
#' @details To avoid issues with the boundary of the parameter support (on the regular scale), the optimization is conducted with a penalty term
#' ```
#' penalty(psi) = log(1+exp(psi[2])) + log(1+exp(-psi[2])) - psi[3].
#' ```
#'
#' @references Savin, T. and Doyle, P.S. "Static and dynamic errors in particle tracking microrheology." Biophysical Journal 88.1 (2005): 623-638.
#'
#' @seealso [fsd_model], the class definition for the fSD model.
#' @example examples/fsd_sim.R
#' @example examples/fsd_fit.R
#'
#' @export
fsd_fit <- function(Xt, dt, drift = c("linear", "none", "quadratic"),
                    vcov = TRUE, ad_only = TRUE) {
  drift <- match.arg(drift)
  model <- fsd_model$new(Xt = Xt, dt = dt, drift = drift)
  # penalty function
  penalty <- function(psi) log1pe(psi[2]) + log1pe(-psi[2]) - psi[3]
  # calculate MLE of penalized profile likelihood
  psi0 <- rep(0, length(model$phi_names))
  opt <- optim(par = psi0,
               fn = function(psi) model$nlp(psi) + penalty(psi))
  if(opt$convergence != 0) warning("`optim()` did not converge.")
  # combine with nuisance parameters
  omega_hat <- model$get_omega(opt$par)
  if(vcov) {
    fi <- model$fisher(omega_hat) # fisher information for likelihood
    fp <- numDeriv::hessian(x = opt$par, func = penalty) # for penalty
    n_phi <- length(model$phi_names)
    fi[1:n_phi,1:n_phi] <- fi[1:n_phi,1:n_phi] + fp
    # convert to (named) variance
    var_hat <- model$get_vcov(fi)
    out <- list(coef = omega_hat, vcov = var_hat)
  } else out <- omega_hat
  if(ad_only) out <- to_ad(out, model, vcov)
  out
  ## model <- fsd_model()
  ## ans <- csi_fit(model, dX, dt, Tz, vcov)
  ## ans
}

# #' @param tau Scalar between 0 and 1, indicating the percentage of time for which the camera shutter is open in the dynamic error model.  Estimated if missing. See Details.
# #' @param sigma2 Magnitude of static error.  Estimated if missing. See Details.

## '
## ' The Savin & Doyle's localization error model with fBM process (fsd model) has the following form:
## ' \deqn{
## ' Y_n = 1/\tau \int_0^\tau X_{n \Delta t +s} ds + \sigma e_{n}
## ' }{
## ' Y[n] = 1/\tau integral_0^\tau X_{n  \Delta t +s} ds + \sigma e[n]
## ' }
## ' where \eqn{X_n} is an `q = 1,2` dimensional fBM model (See [fbm_fit()]).
## '
## ' The integral part measures the average displacement of particle trajectories during \eqn{[n \Delta t, n\Delta t + \tau]} and this term is called dynamic error.
## ' \eqn{e_n} is a standard white noise, and this term is called static error.
## '
## ' When put `tau` = 0, this model becomes fractional static error model.
## ' When put `sigma2` = 0, this model becomes fractional dynamic error model.
## ' When put `(tau, sigma2)` = 0, this model becomes fractional Brownian model.
## '
## ' Optimization is done by [optimize()]. It works better when parameters are re-parametrized into unrestricted form (See [fsd_model()]).
## '


#--- quick test ----------------------------------------------------------------

## fsd_model2 <- function() {
##   # parameter name
##   fsd_theta <- c("alpha", "tau", "sigma2")
##   # ACF function
##   .fsd_acf <- function(theta, dt, N) {
##     fsd_acf(theta[1], theta[2], theta[3], dt, N)
##   }
##   # transformation, from original scale (theta) to unrestricted scale (gamma)
##   fsd_trans <- function(theta) {
##     gamma <- theta
##     # alpha
##     gamma[1] <- logit(theta[1], min = 0, max = 2)
##     # tau
##     gamma[2] <- logit(theta[2], min = 0, max = 1)
##     # sigma2
##     gamma[3] <- log(theta[3])
##     gamma
##   }
##   # inverse transformation, from unrestricted scale (gamma) to original scale (theta)
##   fsd_itrans <- function(gamma) {
##     theta <- gamma
##     # alpha
##     theta[1] <- ilogit(gamma[1], min = 0, max = 2)
##     # tau
##     theta[2] <- ilogit(gamma[2], min = 0, max = 1)
##     # sigma2
##     theta[3] <- exp(gamma[3])
##     theta
##   }
##   # penalty function on transformed scale
##   fsd_penalty <- function(gamma) {
##     log1pe(gamma[2]) + log1pe(-gamma[2]) - gamma[3]
##   }
##   # plug-in values
##   # fsd_plugin <- names_plugin <- NULL
##   # if(!missing(tau)) {
##   #   fsd_plugin <- c(fsd_plugin, tau)
##   #   names_plugin <- c(names_plugin, "tau")
##   # }
##   # if(!missing(sigma2)) {
##   #   fsd_plugin <- c(fsd_plugin, sigma2)
##   #   names_plugin <- c(names_plugin, "sigma2")
##   # }
##   # names(fsd_plugin) <- names_plugin
##   model <- list(theta_names = fsd_theta,
##                 acf = .fsd_acf,
##                 theta_trans = fsd_trans,
##                 theta_itrans = fsd_itrans,
##                 penalty = fsd_penalty)
##   class(model) <- "csi_model"
##   model
## }

## fsd_fit2 <- function(dX, dt, Tz, vcov = TRUE) {
##   model <- fsd_model2()
##   ans <- csi_fit(model, dX, dt, Tz, vcov)
##   ans
## }
