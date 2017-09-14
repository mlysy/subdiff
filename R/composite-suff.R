#' nested function for sufficient statistics of composite likelihood
#' relationship between the nrow of Yt and acf: length(acf) = floor(Yt/ds) - 1
#' @export
composite.suff <- function(Y, X, acf, ds) {
  n <- nrow(Y)
  q <- ncol(Y)
  n.ds <- floor(n / ds) - 1
  
  if(length(X) == 1) {
    X <- matrix(X, n, 1)
  }
  
  noBeta <- all(X == 0)
  
  # variance type
  if(!missing(acf)) {
    if(ncol(acf) != n.ds) {
      stop("Given acf is incompatible with nrow(Y) and npred.")
    }
    if(is.vector(acf)) {
      Tz <- Toeplitz(acf = acf)
    } else if(class(acf) == "Toeplitz") {
      var.type <- "Toeplitz"
      Tz <- acf
      acf <- Tz$getAcf()
    }
  } else {
    stop("Acf missing without default value.")
  }
  
  S <- matrix(0, q, q)
  
  if(!noBeta){
    XvX <- matrix(0, q, q)
    XvY <- matrix(0, q, q)
    
    for(ii in 1:ds) {
      Yt <- downSample(Y, ds, ii)
      Yi <- apply(Yt, 2, diff)
      Xt <- downSample(X, ds, ii)
      Xi <- apply(Xt, 2, diff)
      vX <- solve(Tz, Xi)
      XvX <- XvX + crossprod(Xi, vX)
      XvY <- XvY + crossprod(Yi, vX)
    }
    XvY <- t(XvY)
    Betahat <- solve(XvX, XvY)
  } else {
    Betahat <- matrix(0, 1, q)
  }
  
  for(ii in 1:ds) {
    Yt <- downSample(Y, ds, ii)
    Yi <- apply(Yt, 2, diff)
    Xt <- downSample(X, ds, ii)
    Xi <- apply(Xt, 2, diff)
    Yi <- Yi - Xi %*% Betahat
    S <- S + crossprod(Yi, solve(Tz, Yi))
  }
  ldV <- determinant(Tz)
  
  list(Betahat = Betahat, S = S, ldV = ldV, n.ds = n.ds)
}