// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ModePoly
NumericVector ModePoly(NumericVector roots, double nIter, double tol);
RcppExport SEXP subdiff_ModePoly(SEXP rootsSEXP, SEXP nIterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type roots(rootsSEXP);
    Rcpp::traits::input_parameter< double >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(ModePoly(roots, nIter, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"subdiff_ModePoly", (DL_FUNC) &subdiff_ModePoly, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_subdiff(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
