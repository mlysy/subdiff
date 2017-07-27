//#include <math.h>
//#include <matrix.h>
//#include <mex.h>
#define MACHINE_EPSILON 2.2204460492503131e-16
#include <Rcpp.h>
using namespace Rcpp;
//#ifdef __GNUC__
//#define abs fabs
//#endif

const double seps = sqrt(MACHINE_EPSILON);
const double phi = 2.0/(1.0+sqrt(5.0));

double modePoly(int nIter, double tol, int nRoots, double *roots,
		double lowLim, double hiLim) {
  double tol2, xf, xm;
  double ll, lh, lm1, lm2, Qm1, Qm2;
  int ii, jj;
  // first evaluation
  ll = lowLim;
  lh = hiLim;
  lm1 = ll + (1-phi)*(lh-ll);
  lm2 = ll + phi*(lh-ll);
  Qm1 = 0.0;
  Qm2 = 0.0;
  for(ii = 0; ii < nRoots; ii++) {
    Qm1 -= log(fabs(lm1-roots[ii]));
    Qm2 -= log(fabs(lm2-roots[ii]));
  }
  // subsequent evaluations
  for(jj = 0; jj < nIter; jj++) {
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
        for(ii = 0; ii < nRoots; ii++) {
          Qm1 -= log(fabs(lm1-roots[ii]));
        }
      } else {
        ll = lm1;
        lm1 = lm2;
        Qm1 = Qm2;
        lm2 = ll + phi*(lh-ll);
        Qm2 = 0.0;
        for(ii = 0; ii < nRoots; ii++) {
          Qm2 -= log(fabs(lm2-roots[ii]));
        }
      }
    } else break;
  }
  return(0.5*(ll+lh));
}

//[[Rcpp::export]]
NumericVector ModePoly(NumericVector roots, double nIter, double tol) {
  // local variables
  int nRoots = roots.length();
  // output variables
  NumericVector modes(nRoots-1);
  for(int ii=0; ii < nRoots-1; ii++) {
    modes[ii] = modePoly(nIter, tol, nRoots, REAL(roots),
			 roots[ii], roots[ii+1]);
  }
  return modes;
}

// void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// {
//   // input variables
//   double *roots;
//   double tol;
//   int nIter;
//   // local variables
//   int nRoots, ii;
//   // output variables
//   double *modes;

//   // get inputs
//   roots = mxGetPr(prhs[0]);
//   nRoots = mxGetNumberOfElements(prhs[0]);
//   nIter = (int) *mxGetPr(prhs[1]);
//   tol = (double) *mxGetPr(prhs[2]);

//   // create outputs
//   plhs[0] = mxCreateDoubleMatrix(1, nRoots-1, mxREAL);
//   modes = mxGetPr(plhs[0]);

//   // evaluation
//   for(ii = 0; ii < nRoots-1; ii++) {
//     modes[ii] = modePoly(nIter, tol, nRoots, roots, roots[ii], roots[ii+1]);
//   }

//   return;
// }
