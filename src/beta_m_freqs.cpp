#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston.h"
#include "statistics_class.h"
#include "allelecounter.h"
using namespace Rcpp;
using namespace RcppParallel;

class Betam_freq : public Stats {
  public:
  Betam_freq(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup, IntegerVector Ind_group)
  : Stats(pA, which_snps, SNPgroup, Ind_group) {}

  void compute_stats() {
    if(nb_snps == 0 || nb_snp_groups == 0) {
      return;
    }
    // comptages alléliques
    allelecounter X(&data[0], ncol, true_ncol, nb_snps, nb_ind_groups, ind_group);
    parallelReduce(0, nb_snps, X);
    // std::cout << "parallelReduce ok\n";
 
    // calcul de beta_m avec freqeuences alléliques restimées à chaque permutation
    std::vector<double> U(nb_snps);
    for(size_t i = 0; i < nb_snps; i++) { // boucle sur les SNP
      // frequence allelique dans le groupe 1 (g=0)
      double n = (double) X.R[2*i*nb_ind_groups];
      double m = (double) X.R[2*i*nb_ind_groups + 1];
      double pi = (n+1)/(n+m+2); // comme dans WSS, estimateur "bayésien"
      for(size_t g = 0; g < nb_ind_groups; g++) { // sur les groupes d'invidus
        double n = (double) X.R[2*(i*nb_ind_groups + g)];   // = nb d'alleles 1 dans groupe g+1
        double m = (double) X.R[2*(i*nb_ind_groups + g)+1]; // = nb d'alleles 0 dans groupe g+1
        U[i]    += 0.5*(n*(n-1)/pi + m*(m-1)/(1-pi)) ; 
      }
    }

    // calcule la somme des stats beta_m par région génomique
    for(int i = 0; i < nb_snp_groups; i++) stats[i] = 0;
    for(int i = 0; i < nb_snps; i++) {
      stats[ snp_group[i] - 1 ] +=  U[i]/((double) ncol);
    }
  }

};

//[[Rcpp::export]]
List beta_m_freq(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector region, IntegerVector group, int A_target, int B_max) {
  Betam_freq B(p_A, which_snps, region, group);
  if(B_max > 0) {
    return B.permute_stats(A_target,B_max);
  } else {
    B.compute_stats();
    List L;
    L["statistic"] = B.stats;
    return L;
  }
}

RcppExport SEXP oz_beta_m_freq(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP regionSEXP, SEXP groupSEXP, SEXP A_targetSEXP, SEXP B_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type A_target(A_targetSEXP);
    Rcpp::traits::input_parameter< int >::type B_max(B_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_m_freq(p_A, which_snps, region, group, A_target, B_max));
    return rcpp_result_gen;
END_RCPP
}

