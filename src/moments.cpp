#include <RcppEigen.h>
using namespace Rcpp;

//[[Rcpp::export]]
List moments(NumericMatrix A, NumericMatrix P){
  Eigen::Map<Eigen::MatrixXd> a( as< Eigen::Map<Eigen::MatrixXd> >(A) );
  Eigen::Map<Eigen::MatrixXd> p( as< Eigen::Map<Eigen::MatrixXd> >(P) );
  Eigen::MatrixXd A1 = a*p;
  Eigen::MatrixXd A2 = A1*A1;
  NumericVector c1(4);
  c1[0] = A1.trace();
  c1[1] = A2.trace();

  c1[2] = (A1*A2).trace();
  c1[3] = (A2*A2).trace();
  // faster(?) computation of t(A^3) and t(A^4)
  // THIS IS NOT FASTER
  /* 
  double t3 = 0, t4 = 0;
  int n = A1.rows();
  int m = A1.cols();
  for(int j = 0; j < m; ++j)
    for(int i = 0; i < n; ++i) {
      double a2 = A2(i,j);
      t3 += A1(j,i)*a2;
      t4 += A2(j,i)*a2;
    }
  c1[2] = t3;
  c1[3] = t4;
  */

  double sigmaQ = sqrt(2*c1[1]);
  double s1 = c1[2]/pow(c1[1], 1.5);
  double s2 = c1[3]/pow(c1[1], 2);
  double beta1 = sqrt(8) * s1;
  double beta2 = 12*s2;
  
  List L = List::create(_["mu"] = c1[0], _["sigma"] = sigmaQ, _["skewness"] = beta1, _["kurtosis"] = beta2);
  
  return L;
}

RcppExport SEXP moments(SEXP ASEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    rcpp_result_gen = Rcpp::wrap(moments(A, P));
    return rcpp_result_gen;
END_RCPP
}
