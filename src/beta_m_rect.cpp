#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston.h"
#include "statistics_class.h"
#include "allelecounter.h"
using namespace Rcpp;
using namespace RcppParallel;

class Betam_rect : public Stats {
  public:
  NumericVector full_p;
  std::vector<double> p;

  Betam_rect(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup, IntegerVector Ind_group, NumericVector p_)
  : Stats(pA, which_snps, SNPgroup, Ind_group), full_p(p_), p(full_nb_snps) {
    for(size_t i = 0; i < full_nb_snps; i++) p[i] = full_p[i];
    update_snps(); // (le constructeur Stats() a appelé celui qui est défini dans la classe Stats) 
    // std::cout << "ok\n";
  }

  void compute_stats() {
    if(nb_snps == 0 || nb_snp_groups == 0) {
      return;
    }
    // comptages alléliques
    allelecounter X(&data[0], ncol, true_ncol, nb_snps, nb_ind_groups, ind_group);
    parallelReduce(0, nb_snps, X);
    // std::cout << "parallelReduce ok\n";
 
    // calcul de beta_m rectifié
    std::vector<double> U_phi(nb_snps);
    // std::vector<double> U_pi(nb_snps), J_pi2(nb_snps), J_pi_phi(nb_snps);
    for(size_t i = 0; i < nb_snps; i++) { // boucle sur les SNP
      double pi = p[i];
      if(pi == 0 || pi == 1) continue; // variant monomorphe : ignorer
      for(size_t g = 0; g < nb_ind_groups; g++) { // sur les groupes d'invidus
        double n = (double) X.R[2*(i*nb_ind_groups + g)];   // = nb d'alleles 1 dans groupe g+1
        double m = (double) X.R[2*(i*nb_ind_groups + g)+1]; // = nb d'alleles 0 dans groupe g+1
        U_phi[i]    += 0.5*(n*(n-1)/pi + m*(m-1)/(1-pi) + (n+m)*(n+m-1))/(n+m) ; // le dernier terme pourrait être omis
        // U_pi[i]     += (n/pi - m/(1-pi));
        // J_pi_phi[i] += 0.5*( m*(m-1)/(1-pi)/(1-pi) - n*(n-1)/pi/pi);
        // J_pi2[i]    += -(n/pi/pi+m/(1-pi)/(1-pi));
      }
    }
    // std::cout << "U ok\n";

    // calcule la somme des stats beta_m par région génomique
    for(int i = 0; i < nb_snp_groups; i++) stats[i] = 0;
    for(int i = 0; i < nb_snps; i++) {
      if(p[i] == 0 || p[i] == 1) continue; // variant monomorphe : ignorer
      // SHOW( U_pi[i] )  // --> 0
      // stats[ snp_group[i] - 1 ] += ( U_phi[i] - U_pi[i]*J_pi_phi[i]/J_pi2[i] )/((double) ncol);
      stats[ snp_group[i] - 1 ] +=  U_phi[i]/((double) ncol);
    }
  }

  void update_snps() {
    //Rcout << "derived update\n";
    nb_snps = 0;
    for(bool b : which_snps)
      if(b) nb_snps++;

    data.resize(nb_snps);
    snp_group.resize(nb_snps);
    p.resize(nb_snps);

    for(size_t i = 0; i < nb_snp_groups; i++) 
      nb_snp_in_group[i] = 0;

    // extraction des données pertinentes...
    size_t k = 0;
    for(size_t i = 0; i < full_nb_snps; i++) {
      if(which_snps[i]) {
        p[k] = full_p[i];
        snp_group[k] = full_snp_group[i];
        data[k++] = full_data[i];
        nb_snp_in_group[ full_snp_group[i] - 1 ]++;
      }
    }
  }

};

//[[Rcpp::export]]
List beta_m_rect(XPtr<matrix4> p_A, LogicalVector which_snps, NumericVector p, IntegerVector region, IntegerVector group, int A_target, int B_max) {
  // Rcout << "ça va commencer\n";
  Betam_rect B(p_A, which_snps, region, group, p);
  if(B_max > 0) {
    return B.permute_stats(A_target,B_max);
  } else {
    B.compute_stats();
    List L;
    L["statistic"] = B.stats;
    return L;
  }
}

RcppExport SEXP oz_beta_m_rect(SEXP p_ASEXP, SEXP which_snpsSEXP, SEXP pSexp, SEXP regionSEXP, SEXP groupSEXP, SEXP A_targetSEXP, SEXP B_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type p_A(p_ASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSexp);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< int >::type A_target(A_targetSEXP);
    Rcpp::traits::input_parameter< int >::type B_max(B_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(beta_m_rect(p_A, which_snps, p, region, group, A_target, B_max));
    return rcpp_result_gen;
END_RCPP
}

