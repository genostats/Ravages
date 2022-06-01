#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "GPG.h"

using namespace Rcpp;


// calcule G'PG pour G la matrice bloc-diagonale des ZW ...
// N est le vecteur des tailles de groupes, il faut pondérer par 1/sqrt(ni nj)...
Eigen::MatrixXd GPG(XPtr<matrix4> pA, LogicalVector which_snps, NumericVector p, NumericVector weights, 
                          NumericMatrix P, IntegerVector N, bool symmetrize) {

  std::vector<const uint8_t *> data;
  std::vector<double> pp;
  std::vector<double> W;

  int n  = pA->ncol; // inds
  int m0 = pA->nrow; // snps (avant extraction)
  int k  = N.size(); // groupes
  if(m0 != which_snps.size() || m0 != weights.size() || m0 != p.size()) 
    stop("In GPG, dimensions mismatch");
  if(P.nrow() != P.ncol() || P.nrow() != n*k)
    stop("In GPG, dimensions mismatch (P)");

  // extraction de la sous matrice de SNP considérée
  for(int i = 0; i < m0; i ++) {
    if(which_snps[i]) {
      data.push_back( pA->data[i] );
      pp.push_back( p[i] );
      W.push_back( weights[i] );
    }
  }
  return GPG2(data, pA->true_ncol, pp, W, P, N, symmetrize);
}

