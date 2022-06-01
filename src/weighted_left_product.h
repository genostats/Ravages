#include <Rcpp.h>
#include "gaston/matrix4.h"

using namespace Rcpp;

#ifndef _rvg_wlp_
#define _rvg_wlp_

NumericMatrix WLP(XPtr<matrix4> pA, NumericVector p, const std::vector<double> & we, NumericMatrix & v);
NumericMatrix WLP(const uint8_t ** A_data, const double *p, size_t A_nrow, size_t A_ncol, size_t A_true_ncol, const std::vector<double> & we, NumericMatrix & v);

NumericMatrix WLP2(const uint8_t ** A_data, const double *p, size_t A_nrow, size_t A_ncol, size_t A_true_ncol, const std::vector<double> & we, NumericMatrix & v);

#endif

