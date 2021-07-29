#' @title inference for dynamic error model with fBM process
#' @details 
#' static error model:
#' Y(t) = int_0^tau [G(t-s) * beta + Z(t-s) * Sigma^{1/2}] ds
#' increment form:
#' dY(t) = dG(t) * beta + X(t, tau) * Sigma^{1/2}
#' where X(t, tau) is the integral of increment of fBM process
#' @param Yt Observation of trajectory
#' @param dt interobservation
#' @param acf Toeplitz space
#' @param ConfInt logical, default to be FALSE
#' @return if `ConfInt` = `TRUE` return list containing estimation and length of error bar, 
#' return estimation only otherwise.
#' @note requires package `LMN` and `SuperGauss`
#' @export
dyn.infer <- function(Yt, dt, acf, ConfInt = FALSE) {
  N <- nrow(Yt)
  if(missing(acf)) {
    acf <- Toeplitz(n = N-1)
  } else {
    if(nrow(acf) != (N-1)) {
      stop("Incorrect dimension of acf")
    } 
  }

  est <- optim(par = c(1, .1), fn = dyn.infer.prof, Yt = Yt, dt = dt, acf = acf, 
               control = list(fnscale = -1))$par
  if(ConfInt) {
    err.bar <- hess.estimator(func = dyn.infer.prof, x = est, Yt = Yt, dt = dt, acf = acf)
    list(est = est, err.bar = err.bar)
  } else {
    list(est = est)
  }
}

# nested function for dyn.infer
dyn.infer.prof <- function(theta, Yt, dt, acf){
  if(theta[2] <= 0 || theta[1] <= 0 || theta[2] >= 2) {
    -Inf
  } else {
    alpha <- theta[1]
    tau <- theta[2]
    N <- nrow(Yt)
    d <- ncol(Yt)
    tseq <- 1:N
    t2seq <- 1:N^2
    Gt <- cbind(tseq, t2seq)
    dY <- apply(Yt, 2, diff)
    dG <- apply(Gt, 2, diff)
    V.acf <- fdyn.acf(alpha, tau, dt, N)
    acf$setAcf(V.acf)
    lmn.prof(Y = dY, X = dG, acf = acf)
  }
}
