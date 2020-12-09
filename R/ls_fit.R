#' Fit least-squares estimate.
#'
#' @template args-dX
#' @template args-dt
#' @param lags Integer vector of lags to use in the fit, such that the timepoints used in the fit are `tau = dt * lags`.
#' @param type Either "standard" for the usual LS estimator, or "improved" for the version of Zhang et al (2018).
#' @param vcov If `TRUE`, calculates the variance matrix of `(alpha, logD)`.
#'
#' @return If `vcov = FALSE`, vector of length 2 with estimates of `(alpha, logD)`.  Otherwise, a list with elements `coef` and `vcov`, where the former is the estimate and the latter is the corresponding variance estimator.
#'
#' @note Uses the **subdiff** convention for `D`.
#'
#' @references Zhang, K., Crizer, K.P.R., Schoenfisch, M.H., Hill, B.D., Didier, G. (2018) "Fluid heterogeneity detection based on the asymptotic distribution of the time-averaged mean squared displacement in single particle tracking experiments". *Journal of Physics A: Mathematical and Theoretical*, 51, pp 445601(44).
#' @export
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
  ## tau <- dt * lags
  # need the standard estimator
  y <- log(msd[lags])
  X <- cbind(log(dt * lags), 1)
  qX <- qr(X)
  theta_stand <- solve(qX, y)
  alpha_stand <- theta_stand[1]
  if(type == "improved" || vcov) {
    # variance matrix
    V <- ls_var(alpha_stand, lags, N)
  }
  if(type == "improved") {
    # account for bias
    bias <- ls_bias(alpha_stand, lags, N)
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
#' @param lags Integer vector of length `nlags`.
#' @param N Scalar.
#' @return Bias vector of length `nlags`.
#' @details Output is proportional to formula (18) in Zhang et al (2018), such that
#' ```
#' ls_bias(alpha, lags, N) = lags/N beta_N(alpha, lags).
#' ```
#' @noRd
ls_bias <- function(alpha, lags, N) {
  bias <- 0
  for(ii in (-N+1):(N-1)) {
    bt <- abs(ii/lags + 1)^alpha
    bt <- bt - 2 * abs(ii/lags)^alpha
    bt <- bt + abs(ii/lags - 1)^alpha
    bias <- bias + (1-abs(ii)/N) * bt^2
  }
  .25 * bias/N
}

#' Log-MSD variance matrix.
#'
#' @param alpha Scalar.
#' @param lags Integer vector of length `nlags`.
#' @param N Scalar.
#' @return Variance matrix of size `nlags x nlags`.
#' @details Returns exactly the term `Upsilon(alpha)` in formula (27) of Zhang et al (2018).  This function is now implemented in `src/ls_var.cpp`.
#' @noRd
## ls_var <- function(alpha, lags, N) {
##   tprod <- 1/sqrt(lags %o% lags)
##   t1ov2 <- sqrt(lags %o% (1/lags))
##   t2ov1 <- 1/t1ov2
##   V <- 0
##   for(ii in (-N+1):(N-1)) {
##     Vt <- abs(ii * tprod + t1ov2)^alpha
##     Vt <- Vt - abs(ii * tprod + t1ov2 - t2ov1)^alpha
##     Vt <- Vt - abs(ii * tprod)^alpha
##     Vt <- Vt + abs(ii * tprod - t2ov1)^alpha
##     V <- V + (1-abs(ii)/N) * Vt^2
##   }
##   .5 * V/N
## }
