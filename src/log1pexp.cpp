// numerically stable log(1+exp(x))
// See Malcher (2012): vignette(package = "Rmpfr", topic = "log1mexp-note")

#include <Rcpp.h>
using namespace Rcpp;

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

//[[Rcpp::export("log1pe")]]
NumericVector log1pexp_wrapper(NumericVector x) {
  int n = x.length();
  NumericVector y(n);
  for(int ii=0; ii<n; ii++) {
    y[ii] = log1pexp(x[ii]);
  }
  return y;
}
