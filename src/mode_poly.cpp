/// @file mode_poly.cpp
/// @author Martin Lysy

#define MACHINE_EPSILON 2.2204460492503131e-16
#include <Rcpp.h>
using namespace Rcpp;

/// Calculate a local mode of a polynomial using the golden secant method.
///
/// The polynomial is specified via its roots, such that `p(x) = \prod_{k=1}^K (x - alpha_k)`.
///
/// @param[in] roots Vector containing the roots of `p(x)`.
/// @param[in] n_root Number of roots.
/// @param[in] lower Lower bound for the local mode.
/// @param[in] upper Upper bound for the local mode.
/// @param[in] n_iter Maximum number of iterations.
/// @param[in] tol Relative tolerance.
///
/// @return The value of the mode.
double mode_poly(double *roots, int n_root, 
		 double lower, double upper,
		 int n_iter, double tol) {
  const double seps = sqrt(MACHINE_EPSILON);
  const double phi = 2.0/(1.0+sqrt(5.0));
  double tol2, xf, xm;
  double ll, lh, lm1, lm2, Qm1, Qm2;
  int ii, jj;
  // first evaluation
  ll = lower;
  lh = upper;
  lm1 = ll + (1-phi)*(lh-ll);
  lm2 = ll + phi*(lh-ll);
  Qm1 = 0.0;
  Qm2 = 0.0;
  for(ii = 0; ii < n_root; ii++) {
    Qm1 -= log(fabs(lm1-roots[ii]));
    Qm2 -= log(fabs(lm2-roots[ii]));
  }
  // subsequent evaluations
  for(jj = 0; jj < n_iter; jj++) {
    xm = 0.5*(ll+lh);
    xf = 0.5*(lm1+lm2);
    tol2 = 2.0*(seps*fabs(xf) + tol/3.0);
    if(fabs(xf-xm) > (tol2 - 0.5*(lh-ll))) {
      if(Qm2 > Qm1) {
        lh = lm2;
        lm2 = lm1;
        Qm2 = Qm1;
        lm1 = ll + (1-phi)*(lh-ll);
        Qm1 = 0.0;
        for(ii = 0; ii < n_root; ii++) {
          Qm1 -= log(fabs(lm1-roots[ii]));
        }
      } else {
        ll = lm1;
        lm1 = lm2;
        Qm1 = Qm2;
        lm2 = ll + phi*(lh-ll);
        Qm2 = 0.0;
        for(ii = 0; ii < n_root; ii++) {
          Qm2 -= log(fabs(lm2-roots[ii]));
        }
      }
    } else break;
  }
  return(0.5*(ll+lh));
}

/// Calculate All local modes of a polynomial specified by its roots.
///
/// @param[in] roots Vector containing the polynomial roots.
/// @param[in] n_iter Maximum number of iterations for the golden secant method.
/// @param[in] tol Relative tolerance.
///
/// @return Vector of modes, of length one less than `length(roots)`.
//[[Rcpp::export]]
NumericVector mode_poly(NumericVector roots, int n_iter, double tol) {
  // local variables
  int n_root = roots.length();
  // output variables
  NumericVector modes(n_root-1);
  for(int ii=0; ii < n_root-1; ii++) {
    modes[ii] = mode_poly(REAL(roots), n_root,
			  roots[ii], roots[ii+1],
			  n_iter, tol);
  }
  return modes;
}
