\name{score.reg.mlogit}
\alias{score.reg.mlogit}
\title{
  Logistic/Multinomial regression on genetic score
}
\description{
  Performs a logistical or a non-ordinal multinomial regression on a genetic score
}
\usage{
  score.reg.mlogit(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, burden.score = c("CAST", "WSS", "Other"), other.score = NULL, maf.threshold = 0.01, reflevel, covariates = NULL, alpha = 0.05, get.OR.value = FALSE)
}
\arguments{
  \item{x}{ A bed matrix}
  \item{group}{ A factor containing the group of each individual, x@ped$pheno by default}
  \item{genomic.region}{ A factor containg the genomic region of each SNP, x@snps$genomic.region by default}
  \item{burden.score}{ Which genetic score to use as the dependant variable in the regression: CAST, WSS or another one given by the user}
  \item{other.score}{ A matrix containing the genetic score for each individual (one per row) and for each genomic region (one per column) if it is different from the CAST or the WSS score}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score}
  \item{reflevel}{ The reference group in the regression}
  \item{covariates}{ A matrix containing the potential covariates to include in the model with one covariate per column and one individual per row}
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
  \item{get.OR.value}{ TRUE/FALSE: whether to return the OR values associated to the regression}
}
\details{
  This function is based on the function \code{get.model.parameters.mlogit} which returns the results from each regression.
  If only two groups of individuals are present, a classical logistic regression is performed. If more than two groups of individuals are present, a non-ordinal multinomial regression is performed, comparing each group of individuals to the reference group indicated by the argument \code{reflevel}.
  The genetic score corresponding to CAST or WSS can be directly computed in the function if \code{burden.score}="CAST" or "WSS". If "Other" is indicated, the user needs to give a matrix containing the genetic score he wants to the argument \code{other.score}.This matrix should contain one row per individual and one column per genomic region.
}
\value{
  A dataframe with one row per genomic region and at least two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  \item{OR}{ The OR value(s) associated to the regression, returned if \code{get.OR.value}=TRUE. If there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR}
  \item{l.upper}{ The upper bound of the confidence interval of each OR}
}
\author{
  Ozvan Bocher and Herve Perdry
}
\examples{
  #Examples of the regression on simulated data
  x <- random.bed.matrix.GRR(file.pop.maf = Kryukov, size = c(1000, 500, 500), baseline = c(0.001, 0.001), GRR.matrix = GRR.del, prop.del = 0.5, prop.pro = 0, same.variant = FALSE, genetic.model = "multiplicative", select.gene = "R1", replicates = 5)
  #Filter to keep rare variants
  x <- filter.rare.variants(x, filter = "controls", maf.threshold = 0.01, min.nb.snps = 5)
  #Regression on 3 groups with the score WSS with OR values asked
  score.reg.mlogit(x=x.filter, group=x.filter@ped$pheno, genomic.region=x.filter@snps$genomic.region, burden.score="WSS", maf.threshold=0.01, reflevel="0", alpha=0.05, get.OR.value=TRUE)
}   
