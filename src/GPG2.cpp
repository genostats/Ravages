#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "GPG.h"
#include "weightedGenotypes.h"

using namespace Rcpp;

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

