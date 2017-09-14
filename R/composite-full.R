#' full likelihood for composite
composite.full <- function(Y, X, Beta, Sigma, acf, ds) {
  n <- nrow(Y)
  q <- ncol(Y)
  n.ds <- floor(n / ds) - 1
  
  if(length(X) == 1) {
    X <- matrix(X, n, 1)
  }

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
  
  ldV <- determinant(Tz)
  
  ll <- 0
  for(ii in 1:ds) {
    Yt <- downSample(Y, ds, ii)
    Yi <- apply(Yt, 2, diff)
    Xt <- downSample(X, ds, ii)
    Xi <- apply(Xt, 2, diff)
    Zi <- Yi - Xi %*% Beta
    IP <- tr(solve(Sigma) * crossprod(Zi, solve(Tz, Zi)))
    ll <- ll + n.ds * ldet(Sigma) + IP + q * ldV + n.ds * q * log(2*pi)
  }
  -0.5 * ll
}

ldet <- function(V) {
  determinant(V, log = TRUE)$mod[1]
}

tr <- function(mat) {
  if(length(mat) == 1) {
    ans <- mat
  } else {
    ans <- sum(diag(mat))
  }
  ans
}
