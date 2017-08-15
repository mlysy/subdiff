# Wrapper function for Mode polynomial
# MODEPOLY given a polynomial with distinct roots lambda_0 < ... < lambda_N,
# find its local modes using the golden section method programmed in c++.
# lambda: vector of roots
# nSteps: number of steps
# tol: tolerance
modePoly <- function(lambda, nsteps = 100, tol = 0){
  ModePoly(lambda, nsteps, tol)
}
