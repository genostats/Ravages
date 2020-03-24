\name{burden.mlogit}
\alias{burden.mlogit}
\title{
  Logistic or multinomial regression on a genetic score
}
\description{
  Performs a logistical or a non-ordinal multinomial regression on a genetic score
}
\usage{
burden.mlogit(x, group = x@ped$pheno, 
              genomic.region = x@snps$genomic.region, 
              burden, maf.threshold = 0.01, ref.level, 
              formula, data, get.OR.value = FALSE, 
              alpha = 0.05)
}
\arguments{
  \item{x}{ A bed matrix, only needed if \code{burden="CAST"} or \code{burden="WSS"}}
  \item{group}{ A factor containing the group of each individual, \code{x@ped$pheno} by default}
  \item{genomic.region}{ A factor containg the genomic region of each SNP, \code{x@snps$genomic.region} by default, only needed if \code{burden}="CAST" or \code{burden}="WSS"}
  \item{burden}{"CAST" or "WSS" to directly compute the CAST or the WSS genetic score; or a matrix with one row per individual and one column per \code{genomic.region} if another genetic score is wanted.}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score}
  \item{ref.level}{ The reference group of individuals for the regression}
  \item{formula}{ An R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
  \item{data}{ A matrix containing the potential covariates to include in the model with one covariate per column and one individual per row}
  \item{get.OR.value}{ TRUE/FALSE: whether to return the OR values associated to the regression}
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
}
\details{
  This function will return results from the regression of the phenotype on the genetic score for each genomic region.

  If only two groups of individuals are present, a classical logistic regression is performed. 
  If more than two groups of individuals are present, a non-ordinal multinomial regression is performed, 
  comparing each group of individuals to the reference group indicated by the argument \code{ref.level}. The choice of the reference group won't affect the p-values, but only the Odds Ratios.

  If another genetic score than CAST or WSS is wanted, a matrix with one row per individual and one column per \code{genomic.region} 
  containing this score should be given to \code{burden}. In this situation, no bed matrix \code{x} is needed.
}
\value{
  A dataframe with one row per genomic region and at least two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  If \code{get.OR.value=TRUE}, additional columns are present:
  \item{OR}{ The OR value(s) associated to the regression. If there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR}
  \item{l.upper}{ The upper bound of the confidence interval of each OR}
}
\seealso{ \code{\link{CAST}}, \code{\link{WSS}}, \code{\link{burden.weighted.matrix}} }
\references{ Bocher O, et al. \emph{Rare variant associationtesting for multicategory phenotype.} Genet.Epidemiol. 2019;43:646–656.}
\examples{
#Importation of data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

#Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

#Group variants within known genes
x <- set.genomic.region(x)

#Filter of rare variants: only non-monomorphic variants with
#a MAF lower than 2.5%
#keeping only genomic regions with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

#run burden test WSS, using the 1000Genome population as "outcome"
#The population is available in the object LCT.pop
burden.mlogit(x1, x1@ped$pop, burden = "WSS", ref.level="CEU", get.OR.value=TRUE)

#Simulation of covariates
#with different probabilities in GBR/CEU/FIN and IBS/TSI
covar <- data.frame(sex=c(sample(0:1, sum(table(x1@ped$pop)[c("CEU", "GBR", "FIN")]),TRUE,c(0.2,0.8)), 
                          sample(0:1, sum(table(x1@ped$pop)[c("TSI", "IBS")]), TRUE, c(0.8,0.2))), 
                     u=runif(nrow(x1)))
              
#with the covariate sex
burden.mlogit(x1, x1@ped$pop, burden = "WSS", ref.level="CEU", get.OR.value=TRUE,
              data=covar, formula = ~ sex)
}   
