#' Fit fBM model with composite-likelihood downsampling.
#'
#' @template args-dX
#' @template args-dT
#' @param type type of downsampling method, "naive" means using basic downsampling method, "comp" means using composite likelihood (see Details).
#' @template args-Tz
#' @template args-var_calc
#' @template args-dots_optim
#' @template ret-cov_vcov
#' @details We assume that time series dX follows MN(Beta * dT, V_fBM(alpha, dT), Sigma), thus downsampled dY_ds follows MN(Beta * ds * dT, V_fBM(alpha, dT*ds), Sigma).
#' @export
fds_fit <- function(dX, dT, ds, var_calc = TRUE, Tz) {
  Xt <- apply(dX, 2, cumsum)
  .ds_naive(Xt, dT, Tz, ds, var_calc)
}

.ds_naive <- function(Xt, dT, Tz, ds, var_calc) {
  # downsampling
  Xt <- .down_sample(Xt, ds)
  dX <- apply(Xt, 2, diff)
  fbm_fit(dX, dT*ds, Tz, var_calc)
}

.down_sample <- function(Yt, ds, pos = 1) {
  if(pos > ds) {
    stop("position index must be smaller than downsample rate.")
  }
  if(ds == 1) {
    Yt
  } else {
    N <- nrow(Yt)
    N.ds <- floor(N/ds)
    Yt.ds <- Yt[seq(from = pos, by = ds, length.out = N.ds), ]
    as.matrix(Yt.ds)
  }
}


## .ds.comp <- function(Xt, dT, Tz, ds, var_calc) {
##   # memory allocation
##   N <- nrow(Xt)
##   q <- ncol(Xt)
##   N_ds <- floor(N/ds) - 1
##   nq <- if(q == 1) 1 else 3
##   ntheta <- 1+q+nq
##   theta_hat <- rep(NA, ntheta)
##   theta_names <- c("gamma", paste0("mu", 1:q), paste0("lambda", 1:nq))
##   if(missing(Tz)) Tz <- Toeplitz(n = N_ds)
##   # profile likelihood on transformed scale
##   ll.prof <- function(theta) {
##     alpha <- itrans_alpha(theta)
##     Tz$setAcf(fbm_acf(alpha, dT*ds, N_ds))
##     composite.prof(Y = Xt, X = dT, acf = Tz, ds = ds)
##   }
##   # likelihood on transformed scale
##   loglik <- function(theta) {
##     alpha <- itrans_alpha(theta[1])
##     mu <- theta[1+1:q]
##     Sigma <- itrans_Sigma(theta[1+q+1:nq]) # default: log(D)
##     Tz$setAcf(fbm_acf(alpha, dT*ds, N_ds))
##     composite.full(Y = Xt, X = dT, Beta = mu, Sigma = Sigma, acf = Tz, ds = ds)
##   }
##   # calculate MLE
##   fit <- optimize(f = ll.prof, interval = c(0,2), maximum = TRUE)
##   theta_hat[1] <- fit$maximum # profiled parameters
##   Tz$setAcf(fbm_acf(itrans_alpha(theta_hat[1]), dT*ds, N_ds))
##   suff <- composite.suff(Y = Xt, X = dT, acf = Tz, ds = ds)
##   theta_hat[1+1:q] <- suff$Beta
##   theta_hat[1+q+1:nq] <- trans_Sigma(suff$S/suff$n.ds/ds)
##   names(theta_hat) <- theta_names
##   ans <- theta_hat # no-copy unless ans is modified
##   if(var_calc) {
##     # variance estimate
##     V_hat <- hessian(loglik, x = theta_hat)
##     V_hat <- solveV(-V_hat)
##     colnames(V_hat) <- theta_names
##     rownames(V_hat) <- theta_names
##     ans <- list(coef = theta_hat, vcov = V_hat)
##   }
##   ans
## }

