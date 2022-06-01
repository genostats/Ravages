/*******************************************
 * LOURDEMENT SIMILAIRE A m4_prod_ms.cpp   *
 *                                         *
 * ce produit n'est pas centré             *
 * mais il est pondéré...                  *
 *                                         *
 *******************************************/
#include "weighted_left_product.h"
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;

struct paraWLP2 : public Worker {
  const uint8_t ** data;
  const double * p; // le vecteur de fréquences allélique (freq A2)
  const size_t nrow;
  const size_t ncol;
  const size_t true_ncol;
  const std::vector<double> we;
  const size_t r; // nb cols résultats
  double * v;

  //output
  double * Av;

  //constructeur
  paraWLP2(const uint8_t ** data, const double * p, size_t nrow, size_t ncol, size_t true_ncol, std::vector<double> we, 
          size_t r, double * v, double * Av)
          : data(data), p(p), nrow(nrow), ncol(ncol), true_ncol(true_ncol), we(we), r(r), v(v), Av(Av) { }


  void operator()(size_t beg, size_t end) {
    double gg[4];
    gg[0] = 0;
    for(size_t i = beg; i < end; i++) {
      gg[1] = we[i];
      gg[2] = 2*we[i];
      gg[3] = 2*p[i]*we[i]; // imputation par le "génotype moyen"
      for(size_t c = 0; c < r; c++) {
        size_t k = c*ncol;
        for(size_t j = 0; j < true_ncol; j++) {
          uint8_t x = data[i][j];
          for(int ss = 0; ss < 4 && (4*j + ss < ncol); ss++) {
            Av[c + i*r] += v[k++]*gg[x&3];
            x >>= 2;
          }
        }
      }
    }
  }
};

// calcule R = v' GW // PAS SA TRANSPOSEE
// avec G = (n x m) génotype 0 1 2 (donné par pA), W = matrice diagonale des poids (donnés par we)
// le résultat a dimensions   v.ncol() (r) x nb_snps (m) 
NumericMatrix WLP2(const uint8_t ** A_data, const double * p, size_t A_nrow, size_t A_ncol, size_t A_true_ncol, 
                  const std::vector<double> & we, NumericMatrix & v) {
  int n = A_ncol; // nb inds
  int m = A_nrow; // nb snps
  // Rcout << "n = " << n << " v.nrow = " << v.nrow() << "\n";
  // Rcout << "m = " << m << " we.size = " << we.size() << "\n";
  if(n != v.nrow() || m != we.size()) stop("Dimensions mismatch (WLP2)");
  int r = v.ncol();

  NumericMatrix R(r,m);
  paraWLP2 X(A_data, p, A_nrow, A_ncol, A_true_ncol, we, r, v.begin(), R.begin());

  parallelFor(0, m, X);
  return R;
}

