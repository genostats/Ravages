#include <Rcpp.h>
#include <iostream>
#include "gaston/matrix4.h"
#include "weighted_left_product.h"
#include <RcppEigen.h>

using namespace Rcpp;

#ifndef _RAVAGES_GPG_
#define _RAVAGES_GPG_

// calcule G'PG pour G la matrice bloc-diagonale des ZW ...
// N est le vecteur des tailles de groupes, il faut pondérer par 1/sqrt(ni nj)...
Eigen::MatrixXd GPG(XPtr<matrix4> pA, LogicalVector which_snps, NumericVector p, NumericVector weights, 
                          NumericMatrix P, IntegerVector N, bool symmetrize);

Eigen::MatrixXd GPG(std::vector<const uint8_t *> & data, int true_ncol, std::vector<double> & pp, std::vector<double> & W, NumericMatrix P,
                    IntegerVector N, bool symmetrize);
Eigen::MatrixXd GPG2(std::vector<const uint8_t *> & data, int true_ncol, std::vector<double> & pp, std::vector<double> & W, NumericMatrix P,
                    IntegerVector N, bool symmetrize);
#endif

