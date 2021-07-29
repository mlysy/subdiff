#' Calculate the least-squares estimate of `(alpha, D)`.
#'
#' @name ls_fit
#' @aliases ls_fit ls_msd_fit
#'
#' @template args-Xt
#' @template args-dt
#' @param lags Integer vector of lags to use in the fit, such that the timepoints used in the fit are `tau = dt * lags`.
#' @param type Either "standard" for the usual LS estimator, or "improved" for the version of Zhang et al (2018).
#' @param vcov If `TRUE`, returns an estimate of the variance matrix of `(alpha, logD)`.
#'
#' @return If `vcov = FALSE`, vector of length 2 with estimates of `(alpha, logD)`.  Otherwise, a list with elements `coef` and `vcov`, where the former is the estimate and the latter is the corresponding variance estimator.
#'
#' @details [ls_fit()] first computes the MSD using `msd_fit(Xt, dt, demean = TRUE)` then passes this on to [ls_msd_fit()].  For finer control over the MSD or if it has been precomputed, one may interact with [ls_msd_fit()] directly.
#'
#' @note Uses the **subdiff** convention for `D`.
#'
#' @references Zhang, K., Crizer, K.P.R., Schoenfisch, M.H., Hill, B.D., Didier, G. (2018) "Fluid heterogeneity detection based on the asymptotic distribution of the time-averaged mean squared displacement in single particle tracking experiments". *Journal of Physics A: Mathematical and Theoretical*, 51, pp 445601(44).
#' @export
ls_fit <- function(Xt, dt, lags,
                   type = c("standard", "improved"), vcov = TRUE) {
  # calculate the MSD
  Xt <- check_Xt(Xt)
  N <- nrow(Xt)
  ndim <- ncol(Xt)
  if(max(lags) > N-2) {
    stop("max(lags) must be smaller than nrow(Xt)-2.")
  }
  msd <- msd_fit(Xt = Xt, nlag = max(lags))
  lags <- sort(lags)
  msd <- msd[lags]
  ls_msd_fit(msd = msd, dt = dt, lags = lags, N = N, ndim = ndim,
             type = type, vcov = vcov)
  ## type <- match.arg(type)
  ## has_Xt <- TRUE # depreciated
  ## msd_out <- msd
  ## if(has_Xt) {
  ##   Xt <- check_Xt(Xt)
  ##   # calculate the MSD
  ##   N <- nrow(Xt)
  ##   ndim <- ncol(Xt)
  ##   if(max(lags) > N-2) {
  ##     stop("max(lags) must be smaller than nrow(Xt)-2.")
  ##   }
  ##   msd <- msd_fit(Xt = Xt, nlag = max(lags))
  ##   lags <- sort(lags)
  ##   msd <- msd[lags]
  ## } else {
  ##   N <- Xdim[1]-1
  ##   ndim <- Xdim[2]
  ##   if(length(lags) != length(msd)) {
  ##     stop("msd and lags must have the same length.")
  ##   }
  ##   msd <- msd[order(lags)]
  ##   lags <- sort(lags)
  ## }
  ## ## tau <- dt * lags
  ## # need the standard estimator
  ## y <- log(msd)
  ## X <- cbind(log(dt * lags), 1)
  ## qX <- qr(X)
  ## theta_stand <- solve(qX, y)
  ## alpha_stand <- theta_stand[1]
  ## if(type == "improved" || vcov) {
  ##   # variance matrix
  ##   V <- (1/ndim) * ls_var(alpha_stand, lags, N)
  ## }
  ## if(type == "improved") {
  ##   # account for bias
  ##   bias <- ls_bias(alpha_stand, lags, N)
  ##   suff <- LMN::lmn_suff(Y = as.matrix(y + bias), X = X, V = V)
  ##   theta_hat <- suff$Bhat[,1]
  ##   if(vcov) {
  ##     V_hat <- chol2inv(chol(suff$T))
  ##   }
  ## } else {
  ##   theta_hat <- theta_stand
  ##   if(vcov) {
  ##     # variance of the standard estimator
  ##     V_hat <- solve(qX, t(solve(qX, V)))
  ##   }
  ## }
  ## # format output
  ## tnames <- c("alpha", "logD")
  ## ans <- setNames(theta_hat, tnames)
  ## if(vcov || msd_out) {
  ##   ans <- list(coef = ans)
  ## }
  ## if(vcov) {
  ##   colnames(V_hat) <- rownames(V_hat) <- tnames
  ##   ans <- c(ans, list(vcov = V_hat))
  ## }
  ## if(msd_out) {
  ##   ans <- c(ans, list(msd = msd))
  ## }
  ## ans
}

#' @rdname ls_fit
#'
#' @param msd Vector of empirical MSD estimates computed by [msd_fit()] at the values of `lags`.
#' @param N Length of the particle trajectory.
#' @param ndim Number of dimensions of the particle trajectory.
#' @export
ls_msd_fit <- function(msd, dt, lags, N, ndim,
                       type = c("standard", "improved"), vcov = TRUE) {
  type <- match.arg(type)
  if(length(lags) != length(msd)) {
    stop("msd and lags must have the same length.")
  }
  msd <- msd[order(lags)]
  lags <- sort(lags)
  # first get the standard estimate of (alpha, D)
  y <- log(msd)
  X <- cbind(log(dt * lags), 1)
  qX <- qr(X)
  theta_stand <- solve(qX, y)
  alpha_stand <- theta_stand[1]
  # now improved estimator and/or standard errors
  if(type == "improved" || vcov) {
    # variance matrix
    V <- (1/ndim) * ls_var(alpha_stand, lags, N)
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
  # correct logD by 2*d factor
  theta_hat[2] <- theta_hat[2] - log(2*ndim)
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
