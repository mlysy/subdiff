#' @title Profile Likelihood for Composite
#' @description 
#' @details 
#' @note We assume that the weights [w1, w2, ..., wK] are equal, thus it is omitted in estimation
#' @param Y Position of process
#' @param X Position of mean
#' @param acf Acf of downsampled increment process
#' @param ds Downsampling rate
#' @param noSigma flag for Sigma matrix
#' @return a list with the following elements:
#' \itemize{
#'    \item \code{Betahat = sum(t(Xi)V^{-1}Xi)^{-1}sum(t(Xi)V^{-1}Yi})
#'    \item \code{S = sum(t(Yi-Xi*Betahat)V^{-1}(Yi-Xi*Betahat))}
#'    \item \code{ldV = log(|V|)}
#'    \item \code{n.ds = nrow(Y.ds) - 1}
#'    \item \code{loglik}
#' }
#' @export
composite.prof <- function(suff, Y, X, acf, ds, noSigma = FALSE) {
  # sufficient statistics
  if(missing(suff)) {
    suff <- composite.suff(Y = Y, X = X, acf = acf, ds = ds)
  }
  n.ds <- suff$n.ds
  S <- suff$S
  ldV <- suff$ldV
  q <- nrow(S)
  if(!noSigma) {
    ll <- n.ds*q*(1 + log(2*pi)) + n.ds*ldet(S/n.ds) + q*ldV
  } else {
    ll <- n.ds*q*log(2*pi) + sum(diag(S)) + q*ldV
  }
  -.5 * ds * ll
}

ldet <- function(V) {
  determinant(V, log = TRUE)$mod[1]
}
