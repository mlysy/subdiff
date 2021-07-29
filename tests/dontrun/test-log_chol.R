require(Rcpp)
sourceCpp("log_chol.cpp")

log_chol_R <- function(V) {
  logC <- chol(V)
  diag(logC) <- log(diag(logC))
  logC[upper.tri(logC, diag = TRUE)]
}

ilog_chol_R <- function(logC) {
  # determine size of matrix
  p <- (-1 + sqrt(1 + 8*length(logC)))/2
  C <- matrix(0, p,p)
  C[upper.tri(C, diag = TRUE)] <- logC
  diag(C) <- exp(diag(C))
  crossprod(C)
}


rpd <- function(n) crossprod(matrix(rnorm(n^2), n, n))

n <- sample(1:10, 1)
Sigma <- rpd(n)
range(log_chol(Sigma)- log_chol_R(Sigma))

n <- 3
Sigma <- rpd(n)
nreps <- 1e4
system.time({
  replicate(nreps, log_chol_R(Sigma))
})
system.time({
  replicate(nreps, log_chol(Sigma))
})

n <- 1
nreps <- 1e4
Sigma <- rpd(n)
lambda <- log_chol_R(Sigma)

system.time({
  replicate(nreps, ilog_chol_R(lambda))
})
