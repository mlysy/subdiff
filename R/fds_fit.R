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
fds_fit <- function(dX, dT, ds, method = c("naive", "comp"), 
                    Tz, var_calc = TRUE) {
  Xt <- apply(dX, 2, cumsum)
  if(missing(method)) method = "naive"
  if(method == "naive") {
    Zt <- dsample(Xt, ds)
    dZ <- diff(Zt)
    ans <- fbm_fit(dZ, dT*ds, Tz, var_calc)
  } else {
    # memory allocation
    N <- nrow(dX)
    qq <- ncol(dX)
    nq <- if(qq == 1) 1 else 3
    NN <- floor(N/ds) - 1
    if(missing(Tz)) Tz <- Toeplitz(n = NN)
    theta_hat <- rep(NA, 1+qq+nq)
    theta_names <- c("alpha", paste0("mu", 1:qq), paste0("lambda", 1:nq))
    
    # profile likelihood on transformed scale
    negll.prof <- function(theta) {
      Tz$setAcf(fbm_acf(itrans_alpha(theta), dT*ds, NN))
      nlp <- rep(NA, ds)
      for(ii in 1:ds) {
        Zt <- dsample(Xt, ds, ii)
        dZ <- diff(Zt)
        suff <- lmn.suff(Y = dZ, X = dT*ds, acf = Tz)
        nlp[ii] <- -lmn.prof(suff)
      }
      sum(nlp)
    }
    
    # MLE
    theta_hat[1] <- optimize(f = negll.prof, interval = c(-5, 5))$minimum
    Tz2 <- Toeplitz(acf = fbm_acf(itrans_alpha(theta_hat[1]), dT, N))
    suff <- lmn.suff(Y = dX, X = dT, acf = Tz2)
    theta_hat[1+1:qq] <- suff$Beta
    theta_hat[1+qq+1:nq] <- trans_Sigma(suff$S/suff$n)
    names(theta_hat) <- theta_names
    ans <- theta_hat
    if(var_calc) {
      # likelihood on transformed scale
      negloglik <- function(theta) {
        mu <- theta[1+1:qq]
        Sigma <- itrans_Sigma(theta[1+qq+1:nq])
        Tz2$setAcf(fbm_acf(itrans_alpha(theta[1]), dT, N))
        suff <- lmn.suff(Y = dX, X = dT, acf = Tz2)
        -lmn.loglik(Beta = t(mu), Sigma = Sigma, suff = suff)
      }
      # variance estimate
      V_hat <- hessian(negloglik, x = theta_hat)
      V_hat <- solveV(V_hat)
      colnames(V_hat) <- rownames(V_hat) <- theta_names
      ans <- list(coef = theta_hat, vcov = V_hat)
    }
  }    
  ans
}

dsample <- function(Yt, ds, pos = 1) {
  if(pos > ds) stop("position index must be no greater than ds rate")
  
  as.matrix(Yt[seq(from = pos, by = ds, len = floor(nrow(Yt)/ds)), ])
}
