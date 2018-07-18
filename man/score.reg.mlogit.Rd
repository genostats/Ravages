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
  \item{group}{ A factor containing the group of each individual}
  \item{genomic.region}{ A factor containg the genomic region of each SNP}
  \item{burden.score}{ Which genetic score to use as the dependant variable in the regression}
  \item{other.score}{ A matrix containing the genetic score for each individual if it is different from the CAST or the WSS score}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score}
  \item{reflevel}{ The reference group in the regression, considered as the controls group by default}
  \item{covariates}{ A matrix containing the potential covariates to include in the model with one covariate per column and one individual per row}
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
  \item{get.OR.value}{ TRUE/FALSE: whether to return the OR values associated to the regression}
}
\details{
  This function is based on the function get.model.parameters.mlogit which returns the results from the regression.
  If only two groups of individuals are present, a classical logistic regression is performed. If more than two groups of individuals are present, a non-ordinal multinomial regression is performed, comparing each group of individuals to the reference group indicated by the argument \code{reflevel}.
  The genetic score corresponding to CAST or WSS can be directly computed in the function if \code{burden.score}="CAST" or "WSS". If "Other" is indicated, the user needs to give a matrix containing the genetic score he wants to the argument \code{other.score}.
}
\value{
  A dataframe with one row per genomic region and two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  \item{OR}{ The OR value(s) associated to the regression. If there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR}
  \item{l.upper}{ The upper bound of the confidence interval of each OR}
}
\author{
  Herve Perdry and Ozvan Bocher
}
\examples{
##To fulfill
}   
