#' @title Sufficient Statistics for Profile Composite Likelihood Inference
#'
#' @param Y Original time series following Matrix Normal MN(X*Beta, V, Sigma), downsampled case
#' would be Y_ds ~ MN(X_ds*Beta, V_ds, Sigma)
#' @param X Linear drift of time series. If X is of length 1, X_ds = rep(X, n)
#' @param acf ACF of columnwise-variance matrix \code{V_ds}, either vector or Toeplitz-object
#' @param ds Downsampling rate, starting at 2
#' @note length of acf should equals floor(Yt/ds) - 1
#' @return a list with the following elements:
#' \itemize{
#'    \item \code{Betahat = sum(t(Xi)V^{-1}Xi)^{-1}sum(t(Xi)V^{-1}Yi})
#'    \item \code{S = sum(t(Yi-Xi*Betahat)V^{-1}(Yi-Xi*Betahat))/n.ds/ds}
#'    \item \code{ldV = log(|V|)}
#'    \item \code{n.ds = floor(nrow(Y)/ds) - 1}
#' }
#' @export
composite.suff <- function(Y, X, acf, ds) {
  n <- nrow(Y)
  q <- ncol(Y)
  n.ds <- floor(n / ds) - 1
  noBeta <- all(X == 0)
  if(length(X) == 1) {
    dX <- matrix(X*ds, n.ds, 1)
    p <- 1
  } else {
    p <- ncol(X)
  }

  # variance type
  if(class(acf) == "Toeplitz") {
    if(ncol(acf) != n.ds) {
      stop("Given acf is incompatible with nrow(Y) and ds.")
    }
  } else if(is.vector(acf)) {
    if(length(acf) != n.ds) {
      stop("Given acf is incompatible with nrow(Y) and ds.")
    }
    acf <- Toeplitz(acf = acf)
  } else {
    stop("Acf must be either vector or Toeplitz-class")
  }

  if(!noBeta){
    XvX <- matrix(0, p, p)
    XvY <- matrix(0, q, p)

    for(ii in 1:ds) {
      Yt <- downSample(Y, ds, ii)
      dY <- apply(Yt, 2, diff)
      if(length(X) != 1) {
        Xt <- downSample(Xt, ds, ii)
        dX <- apply(Xt, 2, diff)
      }
      vX <- solve(acf, dX)
      XvX <- XvX + crossprod(dX, vX)
      XvY <- XvY + crossprod(dY, vX)
    }
    XvY <- t(XvY)
    Betahat <- solve(XvX, XvY)
  }

  S <- matrix(0, q, q)
  for(ii in 1:ds) {
    Yt <- downSample(Y, ds, ii)
    dY <- apply(Yt, 2, diff)
    if(!noBeta) {
      if(length(X) != 1) {
        Xt <- downSample(Xt, ds, ii)
        dX <- apply(Xt, 2, diff)
      }
      dY <- dY - dX %*% Betahat
    }
    S <- S + crossprod(dY, solve(acf, dY))
  }
  ldV <- determinant(acf)

  list(Beta = Betahat, S = S, ldV = ldV, n.ds = n.ds)
}
