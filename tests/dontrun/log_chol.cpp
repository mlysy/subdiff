/// @file log_chol.cpp

#include <Rcpp.h>
using namespace Rcpp;

/// Cholesky decomposition of a symmetric, positive definite matrix.
///
/// Calculates the upper triangular Cholesky factor `U` of a symmetric positive definite `n x n` matrix `A`, such that `A = U'U`.
///
/// @param[out] U Matrix of size `n x n` corresponding to the upper triangular Cholesky factor.  The elements below the diagonal remain untouched.  Stored as a vector with matrix elements in column-major order.
/// @param[in] A Symmetric positive-definite matrix of size `n x n`, stored as a vector in column-major order.  The elements below the diagonal are never accessed.
/// @param[in] n Size of the matrices.
///
/// @note The calculation can be performed in-place, i.e., with `U` and `A` referring to the same location in memory.
void chol_decomp(double *U, const double *A, int n) {
  int ii, jj, kk, colI, colJ;
  double tmpSum, tmpInv;
  for(ii = 0; ii < n; ii++) {
    colI = ii*n;
    tmpSum = 0.0;
    for(kk = 0; kk < ii; kk++) {
      tmpSum += U[colI + kk] * U[colI + kk];
    }
    tmpInv = sqrt(A[colI + ii] - tmpSum);
    U[colI + ii] = tmpInv;
    tmpInv = 1.0/tmpInv;
    for(jj = ii+1; jj < n; jj++) {
      colJ = jj*n;
      tmpSum = 0.0;
      for(kk = 0; kk < ii; kk++) tmpSum += U[colJ + kk] * U[colI + kk];
      U[colJ + ii] = tmpInv * (A[colJ + ii] - tmpSum);
    }
  }
  return;
}

/// Compute the log-Cholesky decomposition of a symmetric positive-definite matrix.
// [[Rcpp::export]]
NumericVector log_chol(NumericMatrix Sigma) {
  int n = Sigma.nrow(); // problem dimensions
  NumericVector lambda(n*(n+1)/2); // output i.e., removing zeros
  double* U = new double[n*n]; // temporary storage for Cholesky
  // cholesky decomposition
  chol_decomp(U, REAL(Sigma), n);
  // now walk through it to concatenate elements
  int ii, jj, kk, colI;
  kk = 0;
  for(ii = 0; ii<n; ii++) { // column counter
    colI = ii*n;
    for(jj = 0; jj<=ii; jj++) { // row counter
      lambda[kk] = U[colI + jj];
      if(ii == jj) lambda[kk] = log(lambda[kk]);
      kk++;
    }
  }
  return lambda;
}
