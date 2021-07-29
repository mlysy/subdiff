#include <Rcpp.h>
using namespace Rcpp;

// nested function
// convert MA(1) process into its eps with given rho
//[[Rcpp::export("ma1_resid")]]
NumericMatrix ma1Resid(NumericMatrix Xt, double rho) {
  int N = Xt.nrow();
  int D = Xt.ncol();
  NumericMatrix Yt(N, D);
  double rho1 = 1.0/(1.0-rho);
  for(int jj = 0; jj < D; ++jj){
    Yt(0, jj) = Xt(0, jj) * rho1;
  }
  for(int ii = 1; ii < N; ++ii){
    for(int jj = 0; jj < D; ++jj){
      Yt(ii, jj) = (Xt(ii, jj) - rho * Yt(ii - 1, jj)) * rho1;
    }
  }
  return Yt;
}
