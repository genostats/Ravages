// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston.h"

using namespace Rcpp;
using namespace RcppParallel;



NumericMatrix burden2(matrix4 & A, int n_regions, IntegerVector region, NumericVector w0, NumericVector w1, NumericVector w2) { 
  NumericMatrix B(A.true_ncol*4, n_regions);  
  for(size_t i = 0; i < A.nrow; i++) {
    double ww[4] = {w0(i), w1(i), w2(i), 0};
    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      B(4*j,   region[i]-1) += ww[ (int) d&3 ];
      B(4*j+1, region[i]-1) += ww[ (int) ((d>>2)&3) ];
      B(4*j+2, region[i]-1) += ww[ (int) ((d>>4)&3) ];
      B(4*j+3, region[i]-1) += ww[ (int) ((d>>6)&3) ];
    }
  }
  B = B(Range(0,A.ncol-1), _);
  return B;
}

//[[Rcpp::export]]
NumericMatrix burden2(XPtr<matrix4> p_A, int n_regions, IntegerVector region, NumericVector w0, NumericVector w1, NumericVector w2) {
  return burden2(*p_A, n_regions, region, w0, w1, w2);
}



 
RcppExport SEXP oz_burden2(SEXP p_ASEXP, SEXP n_regionsSEXP, SEXP regionSEXP, SEXP w0SEXP, SEXP w1SEXP, SEXP w2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< int >::type n_regions(n_regionsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w2(w2SEXP);
    rcpp_result_gen = Rcpp::wrap(burden2(p_A, n_regions, region, w0, w1, w2));
    return rcpp_result_gen;
END_RCPP
}


/*
NumericVector burden(matrix4 & A, NumericVector w1, NumericVector w2) { 
  NumericVector B(A.true_ncol*4);  

  for(size_t i = 0; i < A.nrow; i++) {
    double ww[4] = {0, w1(i), w2(i), 0};
    for(size_t j = 0; j < A.true_ncol; j++) {
      uint8_t d = A.data[i][j];
      B[4*j]   += ww[ (int) d&3 ];
      B[4*j+1] += ww[ (int) ((d>>2)&3) ];
      B[4*j+2] += ww[ (int) ((d>>4)&3) ];
      B[4*j+3] += ww[ (int) ((d>>6)&3) ];
    }
  }
  B = B[Range(0,A.ncol-1)];
  return B;
}

NumericVector burden(XPtr<matrix4> p_A, NumericVector w1, NumericVector w2) {
  return burden(*p_A, w1, w2);
}

RcppExport SEXP oz_burden(SEXP p_ASEXP, SEXP w1SEXP, SEXP w2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w1(w1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type w2(w2SEXP);
    rcpp_result_gen = Rcpp::wrap(burden(p_A, w1, w2));
    return rcpp_result_gen;
END_RCPP
}
*/


