/// @file ls_var.cpp

#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;

/// Log-MSD variance matrix.
///
/// @param[in] alpha Scalar.
/// @param[in] tau Vector of length `ntau`.
/// @param[in] N Scalar.
///  
/// @return Variance matrix of size `ntau x ntau`.
/// @details Returns exactly the term `Upsilon(alpha)` in formula (27) of Zhang et al (2018).
// [[Rcpp::export]]
Eigen::MatrixXd ls_var_cpp(double alpha, Eigen::VectorXd tau, int N) {
  int ntau = tau.size();
  // intermediate variables
  MatrixXd tprod = tau * tau.transpose();
  tprod = tprod.cwiseSqrt().cwiseInverse();
  VectorXd itau = tau.cwiseInverse();
  MatrixXd t1ov2 = tau * itau.transpose();
  t1ov2 = t1ov2.cwiseSqrt();
  MatrixXd t2ov1 = t1ov2.transpose();
  MatrixXd Vt = MatrixXd::Zero(ntau, ntau);
  MatrixXd V = MatrixXd::Zero(ntau, ntau);
  // std::cout << "tprod = \n" << tprod << std::endl;
  // std::cout << "t1ov2 = \n" << t1ov2 << std::endl;
  // std::cout << "t2ov1 = \n" << t2ov1 << std::endl;
  // std::cout << "Vt = \n" << Vt << std::endl;
  // MatrixXd Vt2 = V;
  // t1ov2 <- sqrt(tau %o% (1/tau))
  // t2ov1 <- 1/t1ov2
  // V <- 0
  for(int ii=-N+1; ii<=N-1; ii++) {
    // std::cout << "ii = " << ii << std::endl;
    Vt.array() = (ii * tprod + t1ov2).array().abs().pow(alpha);
    // std::cout << "Vt = \n" << Vt << std::endl;
    Vt.array() -= (ii * tprod + t1ov2 - t2ov1).array().abs().pow(alpha);
    // std::cout << "Vt = \n" << Vt << std::endl;
    Vt.array() -= (ii * tprod).array().abs().pow(alpha);
    // std::cout << "Vt = \n" << Vt << std::endl;
    Vt.array() += (ii * tprod - t2ov1).array().abs().pow(alpha);
    // std::cout << "Vt = \n" << Vt << std::endl;
    V.array() += (1.0-abs(1.0 * ii)/N) * Vt.array().square();
    // std::cout << "V = \n" << V << std::endl;
    // Vt <- abs(ii * tprod + t1ov2)^alpha
    // Vt <- Vt - abs(ii * tprod + t1ov2 - t2ov1)^alpha
    // Vt <- Vt - abs(ii * tprod)^alpha
    // Vt <- Vt + abs(ii * tprod - t2ov1)^alpha
    // V <- V + (1-abs(ii)/N) * Vt^2
  }
  V *= .5/N;
  return V;
}
