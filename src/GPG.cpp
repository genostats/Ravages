#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "weighted_left_product.h"
#include <RcppEigen.h>

using namespace Rcpp;


// calcule G'PG pour G la matrice bloc-diagonale des ZW ...
// N est le vecteur des tailles de groupes, il faut pondérer par 1/sqrt(ni nj)...
//[[Rcpp::export]]
NumericMatrix GPG(XPtr<matrix4> pA, LogicalVector which_snps, NumericVector p, NumericVector weights, NumericMatrix P, IntegerVector N) {

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
  int m = 0; // snps (après extraction)
  for(int i = 0; i < m0; i ++) {
    if(which_snps[i]) {
      data.push_back( pA->data[i] );
      pp.push_back( p[i] );
      W.push_back( weights[i] );
      m++;
    }
  }

  // on y va
  NumericMatrix A(k*m, k*m);
  Eigen::Map<Eigen::MatrixXd> AA(as<Eigen::Map<Eigen::MatrixXd>>(A));
  Eigen::Map<Eigen::MatrixXd> PP(as<Eigen::Map<Eigen::MatrixXd>>(P));
  for(int i = 0; i < k; i++) {
    for(int j = i; j < k; j++) {
      Eigen::Block<Eigen::Map<Eigen::MatrixXd>> Pij = PP.block(i*n, j*n, n, n);
      Eigen::Block<Eigen::Map<Eigen::MatrixXd>> Aij = AA.block(i*m, j*m, m, m);
      Eigen::Block<Eigen::Map<Eigen::MatrixXd>> Aji = AA.block(j*m, i*m, m, m);
      NumericMatrix tmp0(n,n); 
      // copie tmp0 = Pij ; on en a besoin car WLP a besoin d'une matrice contigue, pas prévu pour les blocs
      for(int ii = 0; ii < n; ii++) for(int jj = 0; jj < n; jj++) tmp0(jj,ii) = Pij(jj,ii);

      NumericMatrix tmp1 = WLP2(&data[0], &p[0], m, n, pA->true_ncol, W, tmp0);
      NumericMatrix tmp2 = WLP2(&data[0], &p[0], m, n, pA->true_ncol, W, tmp1);
      
      double alpha = 1/sqrt( (double) N[i] * (double) N[j] );
      // il n'y a pas d'opérateur de copie pour assigner Aij = alpha*tmp2
      for(int ii = 0; ii < m; ii++) for(int jj = 0; jj < m; jj++) Aij(jj,ii) = alpha*tmp2(jj,ii);
      Aji = Aij.transpose();
    }
  }
  return A;
}

RcppExport SEXP Ravages_GPG(SEXP pASEXP, SEXP which_snpsSEXP, SEXP pSEXP, SEXP weightsSEXP, SEXP PSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(GPG(pA, which_snps, p, weights, P, N));
    return rcpp_result_gen;
END_RCPP
}
