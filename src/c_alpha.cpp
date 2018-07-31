#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "statistics_class.h"
#include "allelecounter2.h"

using namespace Rcpp;
using namespace RcppParallel;

class Calpha : public Stats {
  public:

  const LogicalVector full_inverse;
  std::vector<bool> inv; // snps à inverser (oblige à redéfinir update_snps)

  IntegerVector full_minor;// alleles mineurs
  std::vector<int> minor;

  IntegerVector full_NNA;// alleles non NA
  std::vector<int> NNA;

  Calpha(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup, IntegerVector ind_group, LogicalVector Inverse)
  : Stats(pA, which_snps, SNPgroup, ind_group), full_inverse(Inverse), full_minor(full_nb_snps), full_NNA(full_nb_snps) {
    if(Inverse.length() != full_nb_snps) 
      stop("Dimensions mismatch");
    // précalculer le nombre de génotypes mineurs et non manquants
    for(size_t i = 0; i < full_nb_snps; i++) {
      for(size_t j = 0; j < true_ncol; j++) {
        uint8_t x = data[i][j];
        for(int ss = 0; ss < 4 && (4*j + ss) < ncol; ss++) {
          if( (x&3) != 3 ) {
            full_NNA[i] += 2;
            if(Inverse[i]) 
              full_minor[i] += 2 -(x&3);
            else
              full_minor[i] += (x&3);
          }
          x >>= 2;
        }
      }
    }
    update_snps(); // le constructeur Stats() a appelé la valeur par défaut
  }

  void update_snps() {
    // Rcout << "derived update\n";
    // nb_snps = sum(which_snps);
    nb_snps = 0;
    for(bool b : which_snps)
      if(b) nb_snps++;

    data.resize(nb_snps);
    snp_group.resize(nb_snps);
    inv.resize(nb_snps); 
    minor.resize(nb_snps);
    NNA.resize(nb_snps);
    for(size_t i = 0; i < nb_snp_groups; i++)
      nb_snp_in_group[i] = 0;

    // extraction des données pertinentes...
    size_t k = 0;
    for(size_t i = 0; i < full_nb_snps; i++) {
      if(which_snps[i]) {
        inv[k] = full_inverse[i];
        NNA[k] = full_NNA[i];
        minor[k] = full_minor[i];
        snp_group[k] = full_snp_group[i];
        data[k++] = full_data[i];
        nb_snp_in_group[ full_snp_group[i] - 1 ]++;
      }
    }
  }

  void compute_stats() {
    if(nb_snps == 0 || nb_snp_groups == 0) {
      return;
    }
    // comptages alléliques
    allelecounter2 X(&data[0], ncol, true_ncol, nb_snps, nb_ind_groups, ind_group, inv);
    parallelReduce(0, nb_snps, X);

    for(int i = 0; i < nb_snp_groups; i++) stats[i] = 0;
    for(size_t i = 0; i < nb_snps; i++) {
      int n = minor[i];
      int grp = snp_group[i] - 1 ;
      for(size_t c = 0; c < nb_ind_groups; c++) {
        int n_c = X.R[ 2*(nb_ind_groups*i+c) ];
        int m_c = X.R[ 2*(nb_ind_groups*i+c) + 1 ];
        double alpha_c = ((double) (n_c+m_c))/((double) NNA[i]);
        stats[ grp ] += (n_c - n*alpha_c)*(n_c - n*alpha_c) - n*alpha_c*(1-alpha_c);
      }
    }

  }


};

//[[Rcpp::export]]
List c_alpha(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector region, IntegerVector group, LogicalVector Inverse, int A_target, int B_max) {

  Calpha B(p_A, which_snps, region, group, Inverse);
  if(B_max > 0) {
    return B.permute_stats(A_target,B_max);
  } else {
    B.compute_stats();
    List L;
    L["statistic"] = B.stats;
    return L;
  }
}

RcppExport SEXP oz_c_alpha(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP regionSEXP, SEXP groupSEXP, SEXP InverseSEXP, SEXP A_targetSEXP, SEXP B_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type Inverse(InverseSEXP);
    Rcpp::traits::input_parameter< int >::type A_target(A_targetSEXP);
    Rcpp::traits::input_parameter< int >::type B_max(B_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(c_alpha(p_A, which_snps, region, group, Inverse, A_target, B_max));
    return rcpp_result_gen;
END_RCPP
}

