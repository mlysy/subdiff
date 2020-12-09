/// @file log1pexp.cpp


#include <Rcpp.h>
using namespace Rcpp;

/// Numerically stable implementation of `log(1 + e^x)`.
///
/// See Malcher (2012): `vignette(package = "Rmpfr", topic = "log1mexp-note")`.
///
/// @param[in] x Real scalar.
/// @return Scalar `y = log(1 + e^x)`.
double log1pexp(double x) {
  double y;
  if(x <= -37.0) {
    y = exp(x);
  } else if(x <= 18.0) {
    y = log1p(exp(x));
  } else if(x <= 33.3) {
    y = x + exp(-x);
  } else {
    y = x;
  }
  return y;
}

/// Numerically stable implementation of `log(1 + e^x)`.
///
/// @param[in] x Input vector of length `n`.
/// @return Vector of length `n` corresponding to `y = log(1 + e^x)`.
///
//[[Rcpp::export]]
NumericVector log1pe(NumericVector x) {
  int n = x.length();
  NumericVector y(n);
  for(int ii=0; ii<n; ii++) {
    y[ii] = log1pexp(x[ii]);
  }
  return y;
}
