#' @title full likelihood for composite
#' @param Y Original time series following Matrix Normal MN(X*Beta, V, Sigma), downsampled case 
#' would be Y_ds ~ MN(X_ds*Beta, V_ds, Sigma)
#' @param X Linear drift of time series. If X is of length 1, X_ds = rep(X, n)
#' @param Beta Parameter for linear drift
#' @param Sigma Parameter for row-wise covariance matrix
#' @param acf ACF of column-wise covariance matrix \code{V_ds}, either vector or Toeplitz-object
#' @param ds Downsampling rate, starting at 2
#' @return Composite likelihood likelihood using profile sufficient statistics
#' @export
composite.full <- function(Y, X, Beta, Sigma, acf, ds) {
  n <- nrow(Y)
  q <- ncol(Y)
  n.ds <- floor(n / ds) - 1
  if(length(X) == 1) {
    dX <- matrix(ds*X, n.ds, 1)
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
  
  ldV <- determinant(acf)
  ldS <- determinant(Sigma, log = TRUE)$mod[1]
  
  ll <- 0
  for(ii in 1:ds) {
    Yt <- downSample(Y, ds, ii)
    dY <- apply(Yt, 2, diff)
    if(length(X) != 1) {
      dX <- downSample(Y, ds, ii)
    }
    dY <- dY - dX %*% Beta
    IP <- sum(diag(solve(Sigma, crossprod(dY, solve(acf, dY)))))
    ll <- ll + n.ds * ldS  + IP + q * ldV
  }
  -0.5 * ll
}
