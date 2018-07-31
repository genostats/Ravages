#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "statistics_class.h"
#include "allelecounter.h"

using namespace Rcpp;
using namespace RcppParallel;

class Betam : public Stats {
  public:

  Betam(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup, IntegerVector ind_group)
  : Stats(pA, which_snps, SNPgroup, ind_group) { }

  void compute_stats() {
    if(nb_snps == 0 || nb_snp_groups == 0) {
      return;
    }
    // comptages alléliques
    allelecounter X(&data[0], ncol, true_ncol, nb_snps, nb_ind_groups, ind_group);
    parallelReduce(0, nb_snps, X);
 
    // calcul de beta_m
    // on calcule directement \sum a_j^2 et \sum a_j ...
    NumericVector R(nb_snps), RR(nb_snps);
    NumericVector S(nb_snps), SS(nb_snps);
    int j = -1;
    for(size_t i = 0; i < nb_ind_groups*nb_snps; i++) {
      if(!(i%nb_ind_groups)) j++;
      R[j]  += X.R[2*i];
      RR[j] += X.R[2*i]*X.R[2*i];
      S[j]  += X.R[2*i+1];
      SS[j] += X.R[2*i+1]*X.R[2*i+1];
    }  
 
    // calcule la somme des stats beta_m par région génomique
    for(int i = 0; i < nb_snp_groups; i++) stats[i] = 0;
    for(int i = 0; i < nb_snps; i++) {
      if(R[i] > 0 && S[i] > 0)
        stats[ snp_group[i] - 1 ] += RR[i]/R[i] + SS[i]/S[i];
    }
  }

};

List beta_m(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector region, IntegerVector group, int A_target, int B_max) {

  Betam B(p_A, which_snps, region, group);
  if(B_max > 0) {
    return B.permute_stats(A_target,B_max);
  } else {
    B.compute_stats();
    List L;
    L["statistic"] = B.stats;
    return L;
  }
}

//[[Rcpp::export]]
List ex_beta_m(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector region, IntegerVector group, IntegerVector g) {

  Betam B(p_A, which_snps, region, group);

  return B.exact_p_value(g);
}



RcppExport SEXP oz_beta_m(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP regionSEXP, SEXP groupSEXP, SEXP A_targetSEXP, SEXP B_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type A_target(A_targetSEXP);
    Rcpp::traits::input_parameter< int >::type B_max(B_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_m(p_A, which_snps, region, group, A_target, B_max));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP oz_ex_beta_m(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP regionSEXP, SEXP groupSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(ex_beta_m(p_A, which_snps, region, group, g));
    return rcpp_result_gen;
END_RCPP
}

