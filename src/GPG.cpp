#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "weighted_left_product.h"
#include <RcppEigen.h>
#include "GPG.h"
#include "weightedGenotypes.h"

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

/*
Eigen::MatrixXd GPG(std::vector<const uint8_t *> & data, int true_ncol, std::vector<double> & pp, std::vector<double> & W, NumericMatrix P, 
                    IntegerVector N, bool symmetrize) {

  int k  = N.size();    // groupes // on ne vérifie plus les dimensions de P 
  int m = data.size();  // snps
  int n = P.ncol() / k; // inds    // confiance totale
//  Rcpp::Rcout << "snps = " << m << "\n";
//  Rcpp::Rcout << "inds = " << n << "\n";
//  Rcpp::Rcout << "gpes = " << k << "\n";
  // on y va
  // NumericMatrix A(k*m, k*m);
  // Eigen::Map<Eigen::MatrixXd> AA(as<Eigen::Map<Eigen::MatrixXd>>(A));
  Eigen::MatrixXd AA(k*m, k*m); 
  Eigen::Map<Eigen::MatrixXd> PP(as<Eigen::Map<Eigen::MatrixXd>>(P));

  for(int i = 0; i < k; i++) {
    for(int j = i; j < k; j++) {
      Eigen::Block<Eigen::Map<Eigen::MatrixXd>> Pij = PP.block(i*n, j*n, n, n);
      Eigen::Block<Eigen::MatrixXd> Aij = AA.block(i*m, j*m, m, m);
      NumericMatrix tmp0(n,n); 
      // copie tmp0 = Pij ; on en a besoin car WLP a besoin d'une matrice contigue, pas prévu pour les blocs
      for(int ii = 0; ii < n; ii++) for(int jj = 0; jj < n; jj++) tmp0(jj,ii) = Pij(jj,ii);

      NumericMatrix tmp1 = WLP2(&data[0], &pp[0], m, n, true_ncol, W, tmp0); // P'G
      NumericMatrix tmp2 = WLP2(&data[0], &pp[0], m, n, true_ncol, W, tmp1); // (P'G)'G = G'PG
      
      double alpha = 1/sqrt( (double) N[i] * (double) N[j] );
      // il n'y a pas d'opérateur de copie pour assigner Aij = alpha*tmp2
      for(int ii = 0; ii < m; ii++) for(int jj = 0; jj < m; jj++) Aij(jj,ii) = alpha*tmp2(jj,ii);

      if(symmetrize) {
        // Eigen::Block<Eigen::Map<Eigen::MatrixXd>> Aji = AA.block(j*m, i*m, m, m);
        Eigen::Block<Eigen::MatrixXd> Aji = AA.block(j*m, i*m, m, m);
        Aji = Aij.transpose();
      }
    }
  }

  return AA;
}
*/

Eigen::MatrixXd GPG2(std::vector<const uint8_t *> & data, int true_ncol, std::vector<double> & pp, std::vector<double> & W, NumericMatrix P, 
                     IntegerVector N, bool symmetrize) {

  int k  = N.size();    // groupes // on ne vérifie plus les dimensions de P 
  int m = data.size();  // snps
  int n = P.ncol() / k; // inds    // confiance totale
//  Rcpp::Rcout << "snps = " << m << "\n";
//  Rcpp::Rcout << "inds = " << n << "\n";
//  Rcpp::Rcout << "gpes = " << k << "\n";
  // on y va
  // NumericMatrix A(k*m, k*m);
  // Eigen::Map<Eigen::MatrixXd> AA(as<Eigen::Map<Eigen::MatrixXd>>(A));
  Eigen::MatrixXd AA(k*m, k*m); 
  Eigen::Map<Eigen::MatrixXd> PP(as<Eigen::Map<Eigen::MatrixXd>>(P));

  Eigen::MatrixXd Gt(weightedGenotypes(&data[0], &pp[0], m, n, true_ncol, W)); // G'

  for(int i = 0; i < k; i++) {
    for(int j = i; j < k; j++) {
      Eigen::Block<Eigen::Map<Eigen::MatrixXd>> Pij = PP.block(i*n, j*n, n, n);
      Eigen::Block<Eigen::MatrixXd> Aij = AA.block(i*m, j*m, m, m);
      // On calcule Aij = G' Pij G / sqrt( ni nj )
      double alpha = 1/sqrt( (double) N[i] * (double) N[j] );
      Aij = alpha * Gt * Pij * Gt.transpose();

      if(symmetrize) {
        Eigen::Block<Eigen::MatrixXd> Aji = AA.block(j*m, i*m, m, m);
        Aji = Aij.transpose();
      }
    }
  }
  return AA;
}

