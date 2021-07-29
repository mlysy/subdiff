#include <Rcpp.h>
using namespace Rcpp;

// nested function
// convert MA(1) process into its eps with given rho
//[[Rcpp::export("ma3_resid")]]
NumericMatrix ma3Resid(NumericMatrix Xt, double rho1, double rho2, double rho3) {
  int N = Xt.nrow();
  int D = Xt.ncol();
  NumericMatrix Yt(N, D);
  double rho = 1.0/(1.0-rho1-rho2-rho3);
  for(int jj = 0; jj < D; ++jj){
    Yt(0, jj) = Xt(0, jj) * rho;
  }
  for(int jj = 0; jj < D; ++jj){
    Yt(1, jj) = (Xt(1, jj) - rho1 * Yt(0, jj)) * rho;
  }
  for(int jj = 0; jj < D; ++jj){
      Yt(2, jj) = (Xt(2, jj) - rho1 * Yt(1, jj) - rho2 * Yt(0, jj)) * rho;
    }
  for(int ii = 3; ii < N; ++ii){
    for(int jj = 0; jj < D; ++jj){
      Yt(ii, jj) = (Xt(ii, jj) - rho1 * Yt(ii-1, jj) - rho2 * Yt(ii-2, jj) - rho3 * Yt(ii-3, jj)) * rho;
    }
  }
  return Yt;
}
