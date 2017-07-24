#include <Rcpp.h>
#include <iostream>
#include "matrix4.h"
#ifndef statclass
#define statclass

class Stats {
  private: 
  Stats();
  public:
  uint8_t ** full_data;
  const int ncol, true_ncol, full_nb_snps;

  const IntegerVector full_snp_group; // facteur pour grouper les SNPs
  const int nb_snp_groups;                // son nombre de niveaux

  LogicalVector which_snps;
  int nb_snps;
  std::vector<uint8_t *> data;
  std::vector<int> snp_group;
  std::vector<int> nb_snp_in_group; // le nombre de snps à TRUE dans which_snp, pour chaque groupe de SNPs

  // output
  NumericVector stats;
  // la fonction qui renvoie un vecteur de stats de longueur nb_snp_groups
  virtual void compute_stats() = 0;

  // la fonction qui permute les phénotypes (les phénotypes et leur 
  // permutations sont définis dans les classes dérivées)
  virtual void permute_pheno() {} ;

  Stats(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup)
    : full_data(pA->data), ncol(pA->ncol), true_ncol(pA->true_ncol), 
    full_nb_snps(pA->nrow), full_snp_group(SNPgroup),
    nb_snp_groups(as<CharacterVector>(SNPgroup.attr("levels")).size()),
    which_snps(which_snps), nb_snp_in_group(nb_snp_groups), stats(nb_snp_groups) {
      if(which_snps.length() != full_nb_snps || SNPgroup.length() != full_nb_snps) {
        stop("Dimensions mismatch\n");
      }
      update_snps();
  }

  void update_snps() {
    nb_snps = sum(which_snps);
    data.resize(nb_snps);
    snp_group.resize(nb_snps);

    for(size_t i = 0; i < nb_snp_groups; i++) 
      nb_snp_in_group[i] = 0;

    // extraction des données pertinentes...
    size_t k = 0;
    for(size_t i = 0; i < full_nb_snps; i++) {
      if(which_snps[i]) {
        snp_group[k] = full_snp_group[i];
        data[k++] = full_data[i];
        nb_snp_in_group[ full_snp_group[i] - 1 ]++;
      }
    }
  }

  List permute_stats(int A_target, int B_max) {
    NumericVector A(nb_snp_groups);
    NumericVector B(nb_snp_groups);
    compute_stats();
    NumericVector Obs = clone(stats);
    // Rcout << "stats = " << stats << "\n";
    for(int b = 0; b < B_max; b++) {
      permute_pheno();
      compute_stats();
      // Rcout << "permutation = " << stats ;
      bool flag = false;
      for(int i = 0; i < nb_snp_groups; i++) {
        if(!nb_snp_in_group[i]) continue;
        B[i]++;
        if(stats[i] >= Obs[i]) {
          A[i]++;
          if(A[i] == A_target) flag = true;
        }
      }
      // Rcout << " A = " << A << " B = " << B << "\n";
      if(!flag) continue;
      // mettre à jour which_snps si nécessaire
      for(int i = 0; i < full_nb_snps; i++) {
        which_snps[i] &= (A[ full_snp_group[i] - 1 ] < A_target); 
      }
      update_snps();
      if(nb_snps == 0) break;
    }
    List L;
    L["stat"] = Obs;
    L["A"] = A;
    L["B"] = B;
    return L;
  }

};


#endif
