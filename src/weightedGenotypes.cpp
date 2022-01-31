#include <iostream>
#include <RcppEigen.h>

// Renvoie (GW)' 
Eigen::MatrixXd weightedGenotypes(const uint8_t ** A_data, const double * p, size_t A_nrow, size_t A_ncol, 
                                 size_t A_true_ncol, const std::vector<double> & W) {
  size_t n = A_ncol; // nb inds
  size_t m = A_nrow; // nb snps
  Eigen::MatrixXd R(m,n);
  
  double gg[4];
  gg[0] = 0;
  for(size_t i = 0; i < m; i++) {
    gg[1] = W[i];
    gg[2] = 2*W[i];
    gg[3] = 2*p[i]*W[i]; // imputation par le "génotype moyen"
    int k = 0;
    for(size_t j = 0; j < A_true_ncol; j++) {
      uint8_t x = A_data[i][j];
      for(int ss = 0; ss < 4 && (4*j + ss < A_ncol); ss++) {
        R(i,k++) = gg[x&3];
        x >>= 2;
      }
    }
  }
  return R;
}

