#' @title inference for static error model with fBM process
#' @details 
#' static error model:
#' Y(t) = G(t) * beta + (Z(t) + a * I(t)) * Sigma^{1/2}
#' increment form:
#' dY(t) = dG(t) * beta + (dZ(t) + a * dI(t)) * Sigma^{1/2}
#' @param Yt Observation of trajectory
#' @param dt interobservation
#' @param acf Toeplitz space
#' @param ConfInt logical, default to be FALSE
#' @return if \code{ConfInt} = \code{TRUE} return list containing estimation and length of error bar, 
#' return estimation only otherwise.
#' @note requires package \code{LMN} and \code{SuperGauss}
#' @export
stat.infer <- function(Yt, dt, acf, ConfInt = FALSE) {
  N <- nrow(Yt)
  if(missing(acf)) {
    acf <- Toeplitz(n = N-1)
  } else {
    if(nrow(acf) != (N-1)) {
      stop("Incorrect dimension of acf")
    } 
  }
  
  est <- optim(par = c(1, .01), fn = stat.infer.prof, Yt = Yt, dt = dt, acf = acf, 
               control = list(fnscale = -1))$par
  if(ConfInt) {
    err.bar <- hess.estimator(func = stat.infer.prof, x = est, Yt = Yt, dt = dt, acf = acf)
    list(est = est, err.bar = err.bar)
  } else {
    list(est = est)
  }
}

# nested function for stat.infer
stat.infer.prof <- function(theta, Yt, dt, acf){
  if(theta[2] < 0 || theta[1] <= 0 || theta[2] >= 2) {
    -Inf
  } else {
    alpha <- theta[1]
    sigma <- sqrt(theta[2])
    N <- nrow(Yt)
    d <- ncol(Yt)
    tseq <- 1:N
    t2seq <- 1:N^2
    Gt <- cbind(tseq, t2seq)
    dY <- apply(Yt, 2, diff)
    dG <- apply(Gt, 2, diff)
    V.acf <- fbm.acf(alpha, dt, N-1) + sigma * c(2, -1, rep(0, N-3))
    acf$setAcf(V.acf)
    lmn.prof(Y = dY, X = dG, acf = acf)
  }
}
