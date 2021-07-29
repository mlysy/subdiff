#include <Rcpp.h>
using namespace Rcpp;

/// Compute the nth power of a polynomial.
///
/// @param[out] bcoef Coefficients of polynomial output.
/// @param[in] acoef Coefficients of polynomial input.
/// @param[in] adeg Degree of input, such that `acoef` is of length `adeg+1`.
/// @param[in] n Power, such that `bcoef` is of length `n*adeg + 1`.
/// @param[in] maxdeg Maximum degree of output polynomial to return.
void poly_pow(double* bcoef, const double* acoef, int adeg, int n, int maxdeg) {
  int bdeg = n * adeg;
  bdeg = (maxdeg < bdeg) ? maxdeg : bdeg;
  int ii,kk;
  double tmp;
  bcoef[0] = std::pow(acoef[0], n);  
  if(bdeg > 0) {
    for(kk=1; kk<=bdeg; kk++) {
      // Rprintf("acoef[%i] = %f\n", kk, acoef[kk]);
      bcoef[kk] = (kk <= adeg) ? n*kk*bcoef[0]*acoef[kk] : 0.0;
      for(ii=1; ii<kk; ii++) {
	tmp = 0.0;
	if(kk-ii <= adeg) tmp += n*bcoef[ii]*acoef[kk-ii];
	if(ii <= adeg) tmp -= acoef[ii]*bcoef[kk-ii];
	bcoef[kk] += (kk-ii) * tmp;
      }
      bcoef[kk] /= kk*acoef[0];
    }
  }
  return;
}

/// Convert AR(p) coefficients to MA(q) approximation.
///
/// The AR(p) model is given by
/// ```
/// X[n] = phi[1] X[n-1] + ... + phi[p] X[n-p] + eps[n].
/// ```
/// The corresponding MA(q) approximation is given by
/// ```
/// X[n] = eps[n] + rho[1] eps[n-1] + ... + rho[q] eps[n-q].
/// ```
/// @param[out] rho Coefficients of MA(q) approximation.
/// @param[in] phi Coefficients of AR(p) model.
/// @param[in] p Number of AR coefficients.
/// @param[in] q Number of MA coefficients.
/// @param[in] bcoef Intermediate storage space.  Must be of length at least `q`.
void ar2ma(double* rho, const double* phi, int p, int q, double* bcoef) {
  int adeg = p-1;
  int bdeg, maxdeg;
  for(int ii=0; ii<q; ii++) rho[ii] = 0.0;
  for(int jj=1; jj<=q; jj++) {
    // Rprintf("jj = %i\n", jj);
    bdeg = adeg*jj;
    maxdeg = q-jj;
    maxdeg = (bdeg < maxdeg) ? bdeg : maxdeg;
    poly_pow(bcoef, phi, adeg, jj, maxdeg);
    // for(int ii=0; ii<=maxdeg; ii++) {
    //   Rprintf("bcoef[%i] = %f\n", ii, bcoef[ii]);
    // }
    for(int ii=0; ii<=maxdeg; ii++) {
      rho[ii+jj-1] += bcoef[ii];
    }
  }
  return;
}

/// Convert ARMA(p,q) coefficients to MA(r) approximation.
///
/// The ARMA(p,q) model is given by
/// ```
/// X[n] = phi[1] X[n-1] + ... + phi[p] X[n-p] + rho[0] eps[n] + ... + rho[q] eps[n-q].
/// ```
/// The corresponding MA(r) approximation is given by
/// ```
/// X[n] = psi[0] eps[n] + psi[1] eps[n-1] + ... + psi[r] eps[n-r].
/// ```
/// @param[out] psi Coefficients of MA(Q) approximation.
/// @param[in] phi Coefficients of autoregressive component of model.
/// @param[in] rho Coefficients of moving-average component of model.
/// @param[in] p Order of autoregressive component of model.
/// @param[in] q Order of moving-average component of model.
/// @param[in] r Order of moving-average approximation.
void arma2ma(double* psi, const double* phi, const double* rho, int p, int q, int r) {
  int max_terms; // maximum number of terms in summation
  int ii, jj;
  for(ii=0; ii<=r; ii++) {
    psi[ii] = ii <= q ? rho[ii] : 0.0;
    max_terms = ii < p ? ii : p;
    for(jj=0; jj<max_terms; jj++) {
      psi[ii] += phi[jj] * psi[ii-jj-1];
    }
  }
  return;
}


// [[Rcpp::export]]
NumericVector poly_pow(NumericVector acoef, int adeg, int n) {
  // int adeg = acoef.length() - 1;
  int bdeg = adeg*n;
  NumericVector bcoef(bdeg + 1);
  poly_pow(REAL(bcoef), REAL(acoef), adeg, n, bdeg);
  return bcoef;
}

// [[Rcpp::export]]
NumericVector ar2ma(NumericVector phi, int q) {
  int p = phi.length();
  NumericVector rho(q);
  // NumericVector bcoef((p-1) * q);
  NumericVector bcoef(q);
  ar2ma(REAL(rho), REAL(phi), p, q, REAL(bcoef));
  return rho;
}

// [[Rcpp::export]]
NumericVector arma2ma(NumericVector phi, NumericVector rho, int q) {
  int p_in = phi.length();
  int q_in = rho.length()-1;
  NumericVector psi(q+1);
  arma2ma(REAL(psi), REAL(phi), REAL(rho), p_in, q_in, q);
  return psi;
}
