#' Fit least-squares estimate.
#'
#' @param dX Matrix of trajectiory increments where each row is an observation and each column is an observed dimension.
#' @param dt Interobservation time.
#' @param lags Integer vector of lags to use in the fit.
#' @param type Either "standard" for the usual LS estimator, or "improved" for the version of Zhang et al (2018).
#' @param vcov If `TRUE`, calculates the variance matrix of `(alpha, logD)`.
#'
#' @return If `vcov = FALSE`, vector of length 2 with estimates of `(alpha, logD)`.  Otherwise, a list with elements `coef` and `vcov`, where the former is the estimate and the latter is the corresponding variance estimator.
#'
#' @note Uses the **subdiff** convention for `D`.
#' @references Zhang, K., Crizer, K.P.R., Schoenfisch, M.H., Hill, B.D., Didier, G. (2018) "Fluid heterogeneity detection based on the asymptotic distribution of the time-averaged mean squared displacement in single particle tracking experiments". *Journal of Physics A: Mathematical and Theoretical*, 51, pp 445601(44).
ls_fit <- function(dX, dt, lags,
                   type = c("standard", "improved"), vcov = TRUE) {
  type <- match.arg(type)
  # calculate the MSD
  N <- nrow(dX)
  lags <- sort(lags)
  if(max(lags) > N-1) {
    stop("max(lags) must be smaller than nrow(dX)-1.")
  }
  msd <- msd_fit(Xt = apply(dX, 2, cumsum), nlag = max(lags))
  tau <- dt * lags
  if(type == "standard" || vcov) {
    # need the standard estimator
    y <- log(msd[lags])
    X <- cbind(1, log(tau))
    qX <- qr(X)
    theta_stand <- solve(qX, y)
    alpha_stand <- theta_stand[2]
  }
  if(type == "improved" || vcov) {
    # variance matrix
    V <- ls_var(alpha_stand, tau, N)
  }
  if(type == "improved") {
    # account for bias
    bias <- ls_bias(alpha_stand, tau, N)
    suff <- LMN::lmn_suff(Y = as.matrix(y + bias), X = X, V = V)
    theta_hat <- suff$Bhat[,1]
    if(vcov) {
      V_hat <- chol2inv(chol(suff$T))
    }
  } else {
    theta_hat <- theta_stand
    if(vcov) {
      # variance of the standard estimator
      V_hat <- solve(qX, t(solve(qX, V)))
    }
  }
  # format output
  tnames <- c("alpha", "logD")
  ans <- setNames(theta_hat, tnames)
  if(vcov) {
    colnames(V_hat) <- rownames(V_hat) <- tnames
    ans <- list(coef = ans, vcov = V_hat)
  }
  ans
}

#--- helper functions ----------------------------------------------------------

#' MSD bias vector.
#'
#' @param alpha Scalar.
#' @param tau Vector of length `ntau`.
#' @param N Scalar.
#' @return Bias vector of length `ntau`.
#' @details Output is proportional to formula (18) in Zhang et al (2018), such that
#' ```
#' ls_bias(alpha, tau, N) = tau/N beta_N(alpha, tau).
#' ```
#' @noRd
ls_bias <- function(alpha, tau, N) {
  bias <- 0
  for(ii in (-N+1):(N-1)) {
    bt <- abs(ii/tau + 1)^alpha
    bt <- bt - 2 * abs(ii/tau)^alpha
    bt <- bt + abs(ii/tau - 1)^alpha
    bias <- bias + (1-abs(ii)/N) * bt^2
  }
  .25 * bias/N
}

#' Log-MSD variance matrix.
#'
#' @param alpha Scalar.
#' @param tau Vector of length `ntau`.
#' @param N Scalar.
#' @return Variance matrix of size `ntau x ntau`.
#' @details Returns exactly the term `Upsilon(alpha)` in formula (27) of Zhang et al (2018).
#' @noRd
ls_var <- function(alpha, tau, N) {
  tprod <- 1/sqrt(tau %o% tau)
  t1ov2 <- sqrt(tau %o% (1/tau))
  t2ov1 <- 1/t1ov2
  V <- 0
  for(ii in (-N+1):(N-1)) {
    Vt <- abs(ii * tprod + t1ov2)^alpha
    Vt <- Vt - abs(ii * tprod + t1ov2 - t2ov1)^alpha
    Vt <- Vt - abs(ii * tprod)^alpha
    Vt <- Vt + abs(ii * tprod - t2ov1)^alpha
    V <- V + (1-abs(ii)/N) * Vt^2
  }
  .5 * V/N
}

#--- scratch -------------------------------------------------------------------

## ls_bias <- function(alpha, tau, N) {
##   iseq <- (-N+1):(N-1)
##   v1 <- 1 - abs(iseq)/N
##   v2 <- abs(iseq/tau+1)^alpha
##   v3 <- abs(iseq/tau)^alpha
##   v3 <- abs(iseq/tau-1)^alpha
##   sum(v1 * (v2 - 2*v3 + v4)^2) / 4*tau
## }

## ls_var <- function(alpha, tau1, tau2, N) {
##   iseq <- (-N+1):(N-1)
##   iseq2 <- iseq/sqrt(tau1*tau2)
##   v1 <- 1 - abs(iseq)/N
##   v2 <- abs(iseq2 + sqrt(tau1/tau2))^alpha
##   v3 <- abs(iseq2 + sqrt(tau1/tau2) - sqrt(tau2/tau1))^alpha
##   v4 <- abs(iseq2)^alpha
##   v5 <- abs(iseq2 - sqrt(tau2/tau1))^alpha
##   sum(v1 * (v2 - v3 - v4 + v5)^2) / 2
## }
