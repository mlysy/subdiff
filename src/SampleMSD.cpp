#include <Rcpp.h>
using namespace Rcpp;

// calculates sample MSD for lag=1,...,nLags.
void sampleMSD(double* MSD, double* X, int N, int nLags) {
  int ii, jj;
  double x2, msd;
  for(ii=0; ii < nLags; ii++) {
    msd = 0.0;
    for(jj=0; jj < N-ii-1; jj++) {
      x2 = X[jj+ii+1] - X[jj];
      msd += x2*x2;
    }
    MSD[ii] += msd/(N-ii-1);
  }
  return;
}

//[[Rcpp::export(".SampleMSD")]]
NumericVector SampleMSD(NumericMatrix X, int nLags) {
  int N = X.nrow();
  int nDims = X.ncol();
  NumericVector msdOut(nLags);
  double *msd = REAL(msdOut);
  double *x = REAL(X);
  for(int ii=0; ii<nDims; ii++) {
    sampleMSD(msd, &x[N*ii], N, nLags);
  }
  return msdOut/(1.0*nDims);
}
