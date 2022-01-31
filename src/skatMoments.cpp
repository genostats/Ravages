#include "GPG.h"

// which_snps permet d'éliminer les SNPs qui n'ont pas la bonne maf etc
// region  : genomic region [FACTEUR]
// group   : groupe des individus [FACTEUR]
// p       : freq alt allele
// Pi      : NullObject$Pi.data 
// P       : NullObject$P1
// weights : SNP weights
// nb_ind_in_group : cf ligne commentées, pour pas le recalcule à chaque fois
// g       : le groupe pour lequel on calcule les moments
// [[Rcpp::export]]
NumericVector skatMoments(XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector region, IntegerVector group,
                 NumericVector p, NumericMatrix Pi, NumericMatrix P, NumericVector weights, 
                 IntegerVector nb_ind_in_group, size_t g) {
  
  // prolegomenes
  size_t nb_snp_groups = as<CharacterVector>(region.attr("levels")).size(); // nb de niveaux
  int m = pA->nrow; // nb SNPS
  if(g >= nb_snp_groups) {
    Rcerr << "g = " << g << "\n";
    stop("[skatMoments] this group does not exist ");
  }
/*
  // pour compter combien d'individus dans chaque groupe
  IntegerVector nb_ind_in_group(nb_ind_groups);
  for(int i : group) nb_ind_in_group[i-1]++;
*/

  std::vector<const uint8_t *> data;
  std::vector<double> pp;
  std::vector<double> W;
  for(int i = 0; i < m; i++) {
    if(which_snps[i] && region[i] == g+1) {
      data.push_back( pA->data[i] );
      pp.push_back( p[i] );
      W.push_back( weights[i] );
    }
  }
  Eigen::MatrixXd A1 = GPG2(data, pA->true_ncol, pp, W, P, nb_ind_in_group, true);

  // calcul des moments à partir de la trace des puissances GPG
  Eigen::MatrixXd A2 = A1*A1;
  NumericVector c1(4); 
  c1[0] = A1.trace();
  c1[1] = A2.trace();
  // Note to self: Eigen does lazy evaluation, there's no possible gain time below
  c1[2] = (A1*A2).trace();
  c1[3] = (A2*A2).trace();

  NumericVector R = NumericVector::create( 
         _["mean"] = c1[0], 
         _["sigma"] = sqrt(2.0 * c1[1]),
         _["skewness"] = 2 * M_SQRT2 * c1[2]/pow(c1[1], 1.5),
         _["kurtosis"] = 12 * c1[3]/(c1[1]*c1[1]) );
 
  return R;
}

RcppExport SEXP rvg_skatMoments(SEXP pASEXP, SEXP which_snpsSEXP, SEXP regionSEXP, SEXP groupSEXP, SEXP pSEXP, SEXP PiSEXP, SEXP PSEXP, SEXP weightsSEXP, SEXP nb_ind_in_groupSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< XPtr<matrix4> >::type pA(pASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type which_snps(which_snpsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type region(regionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group(groupSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Pi(PiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type P(PSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nb_ind_in_group(nb_ind_in_groupSEXP);
    Rcpp::traits::input_parameter< size_t >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(skatMoments(pA, which_snps, region, group, p, Pi, P, weights, nb_ind_in_group, g));
    return rcpp_result_gen;
END_RCPP
}

