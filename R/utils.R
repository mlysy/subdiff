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

#--- get problem dimension -----------------------------------------------------

get_nq <- function(qq) qq*(qq+1)/2

#--- transformations (to be depreciated) ---------------------------------------

trans_alpha <- function(alpha) log(alpha/(2-alpha))
itrans_alpha <- function(gamma) 2/(1+exp(-gamma))
trans_tau <- function(tau) log(tau/(1-tau))
itrans_tau <- function(tau) 1/(1+exp(-tau))
