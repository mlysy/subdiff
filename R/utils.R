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

logit <- function(x) log(x) - log(1-x)
ilogit <- function(x) 1/(1+exp(-x))
