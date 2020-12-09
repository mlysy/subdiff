/// @file msd_empirical.cpp
/// @author Martin Lysy

#include <Rcpp.h>
using namespace Rcpp;

/// Calculate the empirical MSD up to a given number of lags.
///
/// @param[in/out] msd Vector of length `n_lag` containing the computed MSD.  Adds the empirical MSD to whatever `msd` previously contained. 
/// @param[in] Xt Vector containing the trajectory observations (assumed to be equally spaced).
/// @param[in] n_obs Length of trajectory.
/// @param[in] n_lag Number of MSD lags to compute.
void msd_empirical(double* msd, double* Xt, int n_obs, int n_lag) {
  int ii, jj;
  double x2, tmp;
  for(ii=0; ii < n_lag; ii++) {
    tmp = 0.0;
    for(jj=0; jj < n_obs-ii-1; jj++) {
      x2 = Xt[jj+ii+1] - Xt[jj];
      tmp += x2*x2;
    }
    msd[ii] += tmp/(n_obs-ii-1);
  }
  return;
}

/// Calculate the empirical MSD for a multidimensional trajectory.
///
/// @param[in] Xt Matrix where each column is one dimension of the trajectory.
/// @param[in] n_lag Number of MSD lags to compute.
///
/// @return Vector of length `n_lag` containing the *average* MSD per coordinate.
/// 
//[[Rcpp::export]]
NumericVector msd_empirical(NumericMatrix Xt, int n_lag) {
  int n_obs = Xt.nrow();
  int n_dim = Xt.ncol();
  NumericVector msd_out(n_lag);
  double *msd = REAL(msd_out);
  double *xt = REAL(Xt);
  for(int ii=0; ii<n_dim; ii++) {
    msd_empirical(msd, &xt[n_obs*ii], n_obs, n_lag);
  }
  return msd_out/(1.0*n_dim);
}
