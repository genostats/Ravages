#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston.h"
#include "statistics_class.h"

using namespace Rcpp;
using namespace RcppParallel;

struct allelecounter : public Worker {
  // input 
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol; 
  const size_t nrow;
  const size_t nlevels;
  std::vector<int> group; // facteur à nlevels niveaux
  //output
  int * R;

 //constructeur
 allelecounter(uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, const size_t nlevels, std::vector<int> group) 
            : data(data), ncol(ncol), true_ncol(true_ncol), nrow(nrow), nlevels(nlevels), group(group) {
    R = new int[2*nlevels*nrow];
    std::fill(R, R+2*nlevels*nrow, 0);
  }

  //constructeur pour le split
  allelecounter(allelecounter & Q, Split) : data(Q.data), ncol(Q.ncol), true_ncol(Q.true_ncol), nrow(Q.nrow), nlevels(Q.nlevels), group(Q.group) {
    R = new int[2*nlevels*nrow];
    std::fill(R, R+2*nlevels*nrow, 0); 
  }

  // destructeur
  ~allelecounter() {
    delete [] R;
  }

  //worker
  void operator()(size_t beg, size_t end) {
    for(size_t i = beg; i < end; i++) {
      for(size_t j = 0; j < true_ncol; j++) {
        uint8_t x = data[i][j];
        for(int ss = 0; ss < 4 && (4*j + ss) < ncol; ss++) {
          if( (x&3) != 3 ) {
            R[ 2*(i*nlevels + group[4*j+ss]-1) ] += (x&3);
            R[ 2*(i*nlevels + group[4*j+ss]-1) + 1] += 2-(x&3);
          }
          x >>= 2;
        }
      }
    }
  }

  // join
  void join(const allelecounter & Q) {
    std::transform(R, R + 2*nlevels*nrow, Q.R, R, std::plus<int>());
  }

};




class Betam : public Stats {
  public:
  std::vector<int> gr;  // groupe d'individus
  int nlevels;          // nb de groupes

  Betam(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup, IntegerVector group)
  : Stats(pA, which_snps, SNPgroup) {
    nlevels = as<CharacterVector>(group.attr("levels")).size();
    gr.resize(ncol);
    for(size_t i = 0; i < ncol; i++) gr[i] = group[i];
  }

  void permute_pheno() {
    for(int i = ncol - 1; i > 0; i--) {
      int j = (int) std::floor(i*R::runif(0,1));
      int tmp = gr[i];
      gr[i] = gr[j];
      gr[j] = tmp;
    }
  }

  void compute_stats() {
    if(nb_snps == 0 || nb_snp_groups == 0) {
      return;
    }
    // comptages alléliques
    allelecounter X(&data[0], ncol, true_ncol, nb_snps, nlevels, gr);
    parallelReduce(0, nb_snps, X);
 
    // calcul de beta_m
    // on calcule directement \sum a_j^2 et \sum a_j ...
    NumericVector R(nb_snps), RR(nb_snps);
    NumericVector S(nb_snps), SS(nb_snps);
    int j = -1;
    for(size_t i = 0; i < nlevels*nb_snps; i++) {
      if(!(i%nlevels)) j++;
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

//[[Rcpp::export]]
List beta_m(XPtr<matrix4> p_A, LogicalVector which_snps, IntegerVector region, IntegerVector group, int A_target, int B_max) {

  Betam B(p_A, which_snps, region, group);
  if(B_max > 0) {
    return B.permute_stats(A_target,B_max);
  } else {
    B.compute_stats();
    List L;
    L["stat"] = B.stats;
    return L;
  }
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

