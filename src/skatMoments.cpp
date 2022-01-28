#include "GPG.h"

// which_snps permet d'éliminer les SNPs qui n'ont pas la bonne maf etc
// region  : genomic region [FACTEUR]
// group   : groupe des individus [FACTEUR]
// p       : freq alt allele
// Pi      : NullObject$Pi.data 
// P       : NullObject$P1
// weights : SNP weights
// [[Rcpp::export]]
List skatMoments(XPtr<matrix4> pA, LogicalVector which_snps, IntegerVector region, IntegerVector group,
                 NumericVector p, NumericMatrix Pi, NumericMatrix P, NumericVector weights) {
  
  // prolegomenes
  int nb_ind_groups = as<CharacterVector>(group.attr("levels")).size(); // nb de niveaux
  int nb_snp_groups = as<CharacterVector>(region.attr("levels")).size(); // nb de niveaux
  int m = pA->nrow; // nb SNPS
  int n = pA->ncol; // nb inds

  // pour compter combien d'individus dans chaque groupe
  IntegerVector nb_ind_in_group(nb_ind_groups);
  for(int i : group) nb_ind_in_group[i-1]++;

  // Calcul des stats
  // Fabriquer la matrice des Y - Pi
  NumericMatrix Y_Pi(n, nb_ind_groups);
  for(int j = 0; j < nb_ind_groups; j++) { // classes d'individus
    for(int i = 0; i < n; i++) { // individus
      if(group[i] == j+1)
        Y_Pi(i,j) = 1 - Pi(i,j);
      else
        Y_Pi(i,j) = -Pi(i,j);
    }
  }

  // Produit (Y - Pi) * G * W [ dimensions nb_snps x nb_ind_groups : elle est transposée ]
  // extraction matrice "which snps"
  std::vector<const uint8_t *> data;
  std::vector<double> pp;
  std::vector<double> W;
  std::vector<int> region2;
  for(int i = 0; i < m; i++) {
    if(which_snps[i]) {
      data.push_back( pA->data[i] );
      pp.push_back( p[i] );
      W.push_back( weights[i] );
      region2.push_back( region[i] );
    }
  }
  int m1 = data.size(); // nvlle taille
  NumericMatrix Z = WLP(&data[0], &pp[0], m1, n, pA->true_ncol, W, Y_Pi);

  // Vecteur de stats
  NumericVector STATS(nb_snp_groups);
  for(int j = 0; j < nb_ind_groups; j++) {
    for(int i = 0; i < m1; i++) {
      STATS[ region2[i] - 1 ] += Z(i,j)*Z(i,j) / nb_ind_in_group[j];
    }
  }


  // Calcul des moments
  NumericVector MEAN     = no_init(nb_snp_groups);
  NumericVector SIGMA    = no_init(nb_snp_groups);
  NumericVector SKEWNESS = no_init(nb_snp_groups);
  NumericVector KURTOSIS = no_init(nb_snp_groups);

  for(int g = 0; g < nb_snp_groups; g++) {
    // std::vector<const uint8_t *> data;
    // std::vector<double> pp;
    // std::vector<double> W;
    data.clear();
    pp.clear();
    W.clear();
    for(int i = 0; i < m; i++) {
// if(g == 0) Rcout << "(" << (int) which_snps[i] << ":" << region[i] << ")";
      if(which_snps[i] && region[i] == g+1) {
        data.push_back( pA->data[i] );
        pp.push_back( p[i] );
        W.push_back( weights[i] );
      }
    }
    Eigen::MatrixXd A1 = GPG(data, pA->true_ncol, pp, W, P, nb_ind_in_group, true);
//  Rcout << g << " : " << A1.rows() << "x" << A1.cols() << "\n";
    // calcul des moments à partir de la trace des puissances GPG
    Eigen::MatrixXd A2 = A1*A1;
    NumericVector c1(4); 
    c1[0] = A1.trace();
    c1[1] = A2.trace();
    // Note to self: Eigen does lazy evaluation, there's no possible gain time below
    c1[2] = (A1*A2).trace();
    c1[3] = (A2*A2).trace();

    MEAN[g]  = c1[0];
    SIGMA[g] = sqrt(2.0 * c1[1]);
    double s1 = c1[2]/pow(c1[1], 1.5);
    double s2 = c1[3]/(c1[1]*c1[1]);
    SKEWNESS[g] = 2 * M_SQRT2 * s1; // remplace sqrt(8.0) * s1
    KURTOSIS[g] = 12 * s2;
  }

  // On renvoie
  List L;
  L["stat"] = STATS;
  L["mean"] = MEAN;
  L["sigma"] = SIGMA;
  L["skewness"] = SKEWNESS;
  L["kurtosis"] = KURTOSIS;
  return L;
}


RcppExport SEXP rvg_skatMoments(SEXP pASEXP, SEXP which_snpsSEXP, SEXP regionSEXP, SEXP groupSEXP, SEXP pSEXP, SEXP PiSEXP, SEXP PSEXP, SEXP weightsSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(skatMoments(pA, which_snps, region, group, p, Pi, P, weights));
    return rcpp_result_gen;
END_RCPP
}

