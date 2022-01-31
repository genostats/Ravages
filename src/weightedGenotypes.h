using namespace Rcpp;

#ifndef _rvg_wg_
#define _rvg_wg_
// Renvoie (GW)' 
Eigen::MatrixXd weightedGenotypes(const uint8_t ** A_data, const double * p, size_t A_nrow, size_t A_ncol, 
                                 size_t A_true_ncol, const std::vector<double> & W);
#endif

