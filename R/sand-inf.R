#' Downsampling Inference using Composite Likelihood and Sandwich Estimator
#'
#' @param Yt Observation of trajectory
#' @param dt interobservation
#' @param acf Toeplitz space
#' @param rate dowmsapling rate, default to be 1, meaning no downsampling
#' @param ConfInt logical, default to be FALSE
#' @return If \code{ConfInt = TRUE} return list containing estimation and length of error bar, return estimation only otherwise.
#' @note requires package \code{LMN} and \code{SuperGauss}
#' @export
sandwich.infer <- function(Yt, dt, acf, rate = 1, ConfInt = FALSE) {
  N <- nrow(Yt)
  if(missing(acf)) {
    acf <- Toeplitz(n = (floor(N/rate)-1))
  } else {
    if(nrow(acf) != (floor(N/rate)-1)) {
      stop("Incorrect dimension of acf")
    }
  }

  est <- optimize(f = sandwich.infer.prof,
                  interval = c(0, 2), Yt = Yt, dt = dt, acf = acf,
                  rate = rate, control = list(fnscale = -1))$par
  if(ConfInt) {
    err.bar <- sand.estimator(func = sandwich.infer.prof,
                              x = est, Yt = Yt, dt = dt, acf = acf,
                              rate = rate)
    list(est = est, err.bar = err.bar)
  } else {
    list(est = est)
  }
}

# nested function for sandwich.infer
# composite likelihood
sandwich.infer.prof <- function(alpha, Yt, dt, acf, rate) {
  N <- nrow(Yt)
  d <- ncol(Yt)
  N.ds <- floor(N / rate)
  tseq <- 1:N.ds
  t2seq <- (1:N)^2
  Gt <- cbind(tseq, t2seq)
  dG <- apply(Gt, 2, diff)
  V.acf <- fdyn.acf(alpha, tau, dt, N)
  acf$setAcf(V.acf)
  loglik <- 0
  for(ii in 1:rate) {
    Yt.ds <- downSamp(Yt, rate, pos = ii)
    dY.ds <- apply(Yt.ds, 2, diff)
    loglik <- loglik + lmn.prof(Y = dY.ds, X = dG, acf = acf)
  }
  loglik
}

# nested function for sandwich.infer.prof
# downsampling
downSamp <- function(Yt, rate, pos = 1) {
  if(rate == 1){
    Yt
  } else {
    N <- nrow(Yt)
    d <- ncol(Yt)
    N.ds <- floor(N / rate)
    Yt.ds <- matrix(NA, N.ds, d)
    for(ii in 1:N.ds) {
      Yt.ds[ii, ] <- Yt[pos + (ii-1) * rate, ]
    }
    Yt.ds
  }
}
