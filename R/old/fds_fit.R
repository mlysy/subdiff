#' fBM inference using composite downsampling
#' 
#' @param dX one or two-column matrix of trajectory increments.
#' @param dT Interobservation time.
#' @param type type of downsampling method, "naive" means using basic downsampling method, "comp"
#' means using composite downsampling method, "avg" means using average downsampling method
#' @param Tz Optional Toeplitz matrix for intermediate calculations.
#' @param var_calc If \code{TRUE}, also estimate variance matrix.
#' @param ... Additional \code{control} arguments to \code{stats::optim}.
#' @return Vector of coefficients and possibly variance matrix on the transformed scale (see Details).
#' @details We assume that time series dX follows MN(Beta * dT, V_fBM(alpha, dT), Sigma),
#' thus downsampled dY_ds follows MN(Beta * ds * dT, V_fBM(alpha, dT*ds), Sigma).
#' @export
fds_fit <- function(dX, dT, type = "naive", Tz, ds, var_calc = TRUE) {
  Xt <- apply(rbind(0, dX), 2, cumsum)
  switch(type, 
         naive = .ds.naive(Xt, dT, Tz, ds, var_calc = var_calc), 
         comp = .ds.comp(Xt, dT, Tz, ds, var_calc), 
         avg = .ds.avg(Xt, dT, Tz, ds, var_calc))
}

.ds.naive <- function(Xt, dT, Tz, ds, pos = 1, var_calc) {
  # downsampling at position pos
  Xt <- downSample(Xt, ds, pos)
  dX <- apply(Xt, 2, diff)
  # memory allocation and change interobservation time
  N <- nrow(dX)
  q <- ncol(dX)
  dT <- dT*ds
  nq <- if(q == 1) 1 else 3
  ntheta <- 1+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", paste0("mu", 1:q), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N)
  # profile likelihood on transformed scale
  ll.prof <- function(theta) {
    alpha <- itrans_alpha(theta)
    Tz$setAcf(fbm_acf(alpha, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.prof(suff)
  }
  # likelihood on transformed scale
  loglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    mu <- theta[1+1:q]
    Sigma <- itrans_Sigma(theta[1+q+1:nq]) # default: log(D)
    Tz$setAcf(fbm_acf(alpha, dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
    lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
  }
  # calculate MLE
  fit <- optimize(f = ll.prof, interval = c(0, 2), maximum = TRUE)
  theta_hat[1] <- fit$maximum # profiled parameters
  Tz$setAcf(fbm_acf(itrans_alpha(theta_hat[1]), dT, N))
  suff <- lmn.suff(Y = dX, X = dT, acf = Tz)
  theta_hat[1+1:q] <- suff$Beta
  theta_hat[1+q+1:nq] <- trans_Sigma(suff$S/suff$n)
  names(theta_hat) <- theta_names
  ans <- theta_hat # no-copy unless ans is modified
  if(var_calc) {
    # variance estimate
    V_hat <- hessian(loglik, x = theta_hat)
    V_hat <- solveV(-V_hat)
    colnames(V_hat) <- theta_names
    rownames(V_hat) <- theta_names
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}

.ds.avg <- function(Xt, dT, Tz, ds, var_calc) {
  # memory allocation
  N <- nrow(Xt)
  q <- ncol(Xt)
  N_ds <- floor(N/ds) - 1
  nq <- if(q == 1) 1 else 3
  ntheta <- 1+q+nq
  theta_hat <- rep(0, ntheta)
  theta_names <- c("gamma", paste0("mu", 1:q), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N_ds)
  # inference
  if(!var_calc) {
    for(ii in 1:ds) {
      theta_hat <- theta_hat + .ds.naive(Xt, dT, Tz, ds, pos = ii, var_calc)  
    }
    ans <- theta_hat / ds
  } else {
    V_hat <- matrix(0, ntheta, ntheta)
    colnames(V_hat) <- theta_names
    rownames(V_hat) <- theta_names
    for(ii in 1:ds) {
      tmp <- .ds.naive(Xt, dT, Tz, ds, pos = ii, var_calc)
      theta_hat <- theta_hat + tmp$coef
      V_hat <- V_hat + tmp$vcov
    }
    theta_hat <- theta_hat / ds
    V_hat <- V_hat / ds^2
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}

.ds.comp <- function(Xt, dT, Tz, ds, var_calc) {
  # memory allocation
  N <- nrow(Xt)
  q <- ncol(Xt)
  N_ds <- floor(N/ds) - 1
  nq <- if(q == 1) 1 else 3
  ntheta <- 1+q+nq
  theta_hat <- rep(NA, ntheta)
  theta_names <- c("gamma", paste0("mu", 1:q), paste0("lambda", 1:nq))
  if(missing(Tz)) Tz <- Toeplitz(n = N_ds)
  # profile likelihood on transformed scale
  ll.prof <- function(theta) {
    alpha <- itrans_alpha(theta)
    Tz$setAcf(fbm_acf(alpha, dT*ds, N_ds))
    composite.prof(Y = Xt, X = dT, acf = Tz, ds = ds)
  }
  # likelihood on transformed scale
  loglik <- function(theta) {
    alpha <- itrans_alpha(theta[1])
    mu <- theta[1+1:q]
    Sigma <- itrans_Sigma(theta[1+q+1:nq]) # default: log(D)
    Tz$setAcf(fbm_acf(alpha, dT*ds, N_ds))
    composite.full(Y = Xt, X = dT, Beta = mu, Sigma = Sigma, acf = Tz, ds = ds)
  }
  # calculate MLE
  fit <- optimize(f = ll.prof, interval = c(0,2), maximum = TRUE)
  theta_hat[1] <- fit$maximum # profiled parameters
  Tz$setAcf(fbm_acf(itrans_alpha(theta_hat[1]), dT*ds, N_ds))
  suff <- composite.suff(Y = Xt, X = dT, acf = Tz, ds = ds)
  theta_hat[1+1:q] <- suff$Beta
  theta_hat[1+q+1:nq] <- trans_Sigma(suff$S/suff$n.ds/ds)
  names(theta_hat) <- theta_names
  ans <- theta_hat # no-copy unless ans is modified
  if(var_calc) {
    # variance estimate
    V_hat <- hessian(loglik, x = theta_hat)
    V_hat <- solveV(-V_hat)
    colnames(V_hat) <- theta_names
    rownames(V_hat) <- theta_names
    ans <- list(coef = theta_hat, vcov = V_hat)
  }
  ans
}

