#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"
#include "polytomic_haplotype_proba.h"

using namespace Rcpp;

// haplos = matrice d'haplotypes (dimensions u x p) (u haplotypes, p SNPs)
// freq = matrices de fréquences haplotypiques dans les divers groupes
// (dim u x K) (K groupes d'individus)
// size un vecteur de tailles (longueur K)
// reps = nombre de réplicats 
// [[Rcpp::export]]
XPtr<matrix4> rbm_haplos_freqs(NumericMatrix haplos, NumericMatrix freq, NumericVector size, int reps) {
  int u = haplos.nrow(); // nb d'haplotypes
  int p = haplos.ncol(); // nb SNPs
  int K = freq.ncol(); // nb groupes de pop
  if(freq.nrow() != u || size.length() != K)
    stop("Dimensions mismatch");

  int N = sum(size);
  XPtr<matrix4> pA(new matrix4(reps*p, N));

  for(int r = 0; r < reps; r++) { // replicats
    // tirage liste d'haplotypes pour le replicat
    // each elt is a vector of 2*size[k] haplotypes for the size[k] individuals of group k = 0 .. K-1
    List LHaps(K);
    for(int k = 0; k < K; k++) { // groupe de pop k
      NumericVector p = freq(_ , k);  // freqs haplo groupe k;
      LHaps[k] = sample(u, 2*size[k], true, p, false);
    }
    // filling data for current replicate, SNP after SNP
    for(int snp = 0; snp < p; snp++) {
      int ind = 0; // indice de l'individu
      for(int k = 0; k < K; k++) { // groupe de pop k 
        IntegerVector H = LHaps[k];
        for(int i = 0; i < size[k]; i++) {
          pA->set(r*p + snp, ind++, haplos(H[2*i],snp)+haplos(H[2*i+1],snp) );
        }
      }
    }
  }
  return pA;
}

RcppExport SEXP rbm_haplos_freqs(SEXP haplosSEXP, SEXP freqSEXP, SEXP sizeSEXP, SEXP repsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type haplos(haplosSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< int >::type reps(repsSEXP);
    rcpp_result_gen = Rcpp::wrap(rbm_haplos_freqs(haplos, freq, size, reps));
    return rcpp_result_gen;
END_RCPP
}
