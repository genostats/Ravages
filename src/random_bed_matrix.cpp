#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "gaston.h"

using namespace Rcpp;

//[[Rcpp::export]]
XPtr<matrix4> random_bed_matrix(NumericMatrix maf, NumericVector size) {
  int nb_pop = maf.nrow();
  if(size.length() != nb_pop)
    stop("Dimensions mismatch");
  int nrow = maf.ncol(); // nb snp
  int ncol = sum(size);  // nb inds
  XPtr<matrix4> p_A(new matrix4(nrow, ncol));
  for(int i = 0; i < nrow; i++) {
    int k = 0;
    for(int pop = 0; pop < nb_pop; pop++) {
      double q = maf(pop, i);
      double pr1 = (1-q)*(1-q), pr2 = pr1 + 2*(1-q)*q;
      int s = size[pop];
      for(int j = 0; j < s; j++) {
        double r = Rf_runif(0.0, 1.0);
        if(r < pr1) { 
          p_A->set(i,k++,0);
        } else if(r < pr2) {
          p_A->set(i,k++,1);
        } else {
          p_A->set(i,k++,2);
        }
      }
    }
  }
  return p_A;
}

RcppExport SEXP oz_random_bed_matrix(SEXP mafSEXP, SEXP sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type maf(mafSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(random_bed_matrix(maf, size));
    return rcpp_result_gen;
END_RCPP
}

