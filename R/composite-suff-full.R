#' composite full-suff
#' @export
composite.suff.full <- function(Y, X, acf, ds) {
  N <- nrow(Y)
  N.ds <- floor(N/ds) - 1
  suff <- composite.suff(Y, X, acf, ds)
  Beta <- suff$Betahat
  Sigma <- suff$S/N.ds/ds
  composite.full(Y, X, Beta, Sigma, acf, ds)
}