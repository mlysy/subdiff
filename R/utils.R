#--- solve method for variance matrices ----------------------------------------

# optionally computes log determinant as well
solveV <- function(V, x, ldV = FALSE) {
  C <- chol(V)
  if(missing(x)) x <- diag(nrow(V))
  ans <- backsolve(r = C, x = backsolve(r = C, x = x, transpose = TRUE))
  if(ldV) {
    ldV <- 2 * sum(log(diag(C)))
    ans <- list(y = ans, ldV = ldV)
  }
  ans
}

#--- logit and its inverse -----------------------------------------------------

logit <- function(x, min = 0, max = 1) {
  x <- (x - min) / (max - min)
  log(x) - log(1-x)
}
ilogit <- function(x, min = 0, max = 1) {
  1/(1+exp(-x)) * (max - min) + min
}

#--- get problem dimension -----------------------------------------------------

get_nq <- function(qq) {
  if(qq == 1) ans <- 1
  if(qq == 2) ans <- 3
  ans
}
