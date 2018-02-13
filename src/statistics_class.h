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

  LogicalVector which_snps_orig;
  std::vector<bool> which_snps;
  int nb_snps;
  std::vector<uint8_t *> data;
  std::vector<int> snp_group;
  std::vector<int> nb_snp_in_group; // le nombre de snps à TRUE dans which_snp, pour chaque groupe de SNPs
  
  // pour les calculs exacts
  std::vector<int> no_var, some_var;

  // output
  NumericVector stats;
  // la fonction qui renvoie un vecteur de stats de longueur nb_snp_groups
  virtual void compute_stats() = 0;

  // la fonction qui permute les phénotypes (les phénotypes et leur 
  // permutations sont définis dans les classes dérivées)
  virtual void permute_pheno() {} ;

  Stats(const XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector SNPgroup) : 
    full_data(pA->data), 
    ncol(pA->ncol), 
    true_ncol(pA->true_ncol), 
    full_nb_snps(pA->nrow), 
    full_snp_group(SNPgroup),
    nb_snp_groups(as<CharacterVector>(SNPgroup.attr("levels")).size()),
    which_snps_orig(which_snps), 
    which_snps(full_nb_snps),
    nb_snp_in_group(nb_snp_groups), 
    stats(nb_snp_groups)
  {
      if(which_snps_orig.length() != full_nb_snps || SNPgroup.length() != full_nb_snps) {
        stop("Dimensions mismatch\n");
      }
      for(int i = 0; i < full_nb_snps; i++)
        which_snps[i] = which_snps_orig[i];
      update_snps();
  }

  virtual void update_snps() {
    //Rcout << "original update\n";
    // count 'true' in which_snps
    nb_snps = 0;
    for(bool b : which_snps)
      if(b) nb_snps++;

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
    IntegerVector A(nb_snp_groups);
    IntegerVector B(nb_snp_groups);
    IntegerVector C(nb_snp_groups);
    compute_stats();
    NumericVector Obs = clone(stats);
    // Rcout << "stats = " << stats << "\n";
    for(int b = 0; b < B_max; b++) {
      permute_pheno();
      compute_stats();
      // Rcout << "permutation = " << stats ;
      bool flag = false; // ce drapeau se lèvera si A_target est atteint dans un groupe
      for(int i = 0; i < nb_snp_groups; i++) {
        if(!nb_snp_in_group[i]) continue;
        B[i]++;
       
        if(stats[i] >= Obs[i]) {
          A[i]++;
          if(A[i] == A_target) flag = true;
          if(stats[i] == Obs[i]) {
            C[i]++;
          } 
        } 
      }
      // Rcout << " A = " << A << " B = " << B << "\n";
      if(!flag) continue;
      // mettre à jour which_snps si nécessaire (le drapeau est levé !)
      for(int i = 0; i < full_nb_snps; i++) {
        which_snps[i] = which_snps[i] && (A[ full_snp_group[i] - 1 ] < A_target); 
      }
      update_snps();
      if(nb_snps == 0) break;
    }
    List L;
    L["stat"] = Obs;
    L["nb.geq"] = A;
    L["nb.eq"] = C;
    L["nb.perms"] = B;
    NumericVector p(nb_snp_groups);
    for(int j = 0; j < nb_snp_groups; j++) {
      p[j] = ((double) A[j] - 0.5* (double) C[j] + 1.0)/((double) B[j] + 1.0);
    }
    L["p.value"] = p;
    return L;
  }

  /**********************************************
   * des fonctions pour les tests exacts
   **********************************************/

  void keep_one_snp_group(int group) {
    for(size_t i = 0; i < nb_snp_groups; i++) 
      which_snps[i] = which_snps_orig[i] && ( full_snp_group[i] == group );
    update_snps();
  }

  // quels sont les individus qui portent des variants rares ?
  void set_non_zero_inds() {
    // stat Individus
    std::vector<int> stat_inds(16*true_ncol);
    for(uint8_t * da : data) {
      for(size_t j = 0; j < true_ncol; j++) {
        uint8_t d = da[j];
        stat_inds[16*j + ((int) d&3)]++;
        stat_inds[16*j + 4 + (int) ((d>>2)&3)]++;
        stat_inds[16*j + 8 + (int) ((d>>4)&3)]++;
        stat_inds[16*j + 12 + (int) ((d>>6)&3)]++;
      }
    }

    no_var.clear();
    some_var.clear();
    for(size_t j = 0; j < ncol; j++) {
      int N0  = stat_inds[4*j];
      int N1  = stat_inds[4*j+1];
      int N2  = stat_inds[4*j+2];
      int NAs = stat_inds[4*j+3];
      if(N1 == 0 && (N2 == 0 || N0 == 0))
        no_var.push_back(j);
      else
        some_var.push_back(j);
    }
  }
};


#endif
