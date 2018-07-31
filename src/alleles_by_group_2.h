#include <Rcpp.h>
#include <RcppParallel.h>
#include <iostream>
#include <ctime>

#ifndef plic
#define plic
using namespace Rcpp;
using namespace RcppParallel;

struct plic : public Worker {
  // input 
  uint8_t ** data;
  const size_t ncol;
  const size_t true_ncol; 
  const size_t nrow;
  const size_t nlevels;
  std::vector<int> group; // facteur à nlevels niveaux
  std::vector<bool> inverse;
  //output
  int * R;

  //constructeur
  plic(uint8_t ** data, const size_t ncol, const size_t true_ncol, const size_t nrow, const size_t nlevels, std::vector<int> group, std::vector<bool> inverse);
  //constructeur pour le split
  plic(plic & Q, Split);
  // destructeur
  ~plic();
  //worker
  void operator()(size_t beg, size_t end);
  // join
  void join(const plic & Q);
};
#endif
