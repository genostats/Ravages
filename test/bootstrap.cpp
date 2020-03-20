//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

// Générer des Y avec la bonne proba, une ligne de Pi après l'autre
std::pair<NumericMatrix,NumericVector> random_group(NumericMatrix Pi) {
  int n = Pi.nrow(), m = Pi.ncol();
  NumericMatrix Y(n,m);
  NumericVector u( runif(n, 0., 1.) );
  NumericVector cp(n);
  NumericVector nb_inds(m);
  for(int j = 0; j < m; ++j) 
    for(int i = 0; i < n; ++i)
      if( u[i] >= cp[i] && (cp[i] += Pi(i,j)) > u[i] ) {
        Y(i,j) = 1;
        nb_inds(j)++;
      }
  return std::make_pair(Y, nb_inds);
}

// renvoie une matrice de (y - pi^) bootstrapés + le nb de cas par catégorie
// U = Id - U par rapport à bootstrap.r
//[[Rcpp::export]]
List bootstrap1(NumericMatrix U, NumericMatrix Pi) {
  Eigen::Map<Eigen::MatrixXd> u( as<Eigen::Map<Eigen::MatrixXd>>(U) );
  int n = Pi.nrow(), m = Pi.ncol();
  auto rand(random_group(Pi));
  NumericMatrix Y = rand.first;

  NumericMatrix ymp(n,m);
  // Y - pi^_0 où pi^_0  est le vecteur des pi^ estimé sur les données
  // les termes +n pemettent de supprimer la premiere colonne de Y et Pi 
  Eigen::VectorXd ymp0 = Eigen::Map<Eigen::VectorXd>( &Y[0] + n, n*(m-1) ) - Eigen::Map<Eigen::VectorXd>( &Pi[0] + n, n*(m-1) ) ;
  // Rcout << ymp0; 
  // ymp[,-1] vue comme vecteur
  Eigen::Map<Eigen::VectorXd> ymp_( &ymp[0] + n, n*(m-1) );
  ymp_ = u * ymp0;
  // compléter la colonne 0
  for(int j = 1; j < m; ++j) 
    for(int i = 0; i < n; ++i)
       ymp(i,0) -= ymp(i,j);
 
  List L;
  L["ymp"] = ymp;
  L["n"] = rand.second; 
  return L;
}


