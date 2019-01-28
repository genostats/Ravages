#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>
#include "gaston/matrix4.h"

using namespace Rcpp;
using namespace RcppParallel;

NumericMatrix WLP(XPtr<matrix4> pA, const std::vector<double> & we, NumericMatrix & v);
NumericMatrix WLP(const uint8_t ** A_data, size_t A_nrow, size_t A_ncol, size_t A_true_ncol, const std::vector<double> & we, NumericMatrix & v);

