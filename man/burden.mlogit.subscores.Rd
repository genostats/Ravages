\name{burden.mlogit.subscores}
\alias{burden.mlogit.subscores}
\title{
  Logistic or multinomial regression on a multiple genetic scores within a genomic region
}
\description{
  Performs burden tests with subscores in the regression on categorical phenotypes
}
\usage{
burden.mlogit.subscores(x, NullObject, genomic.region = x@snps$genomic.region, 
                        SubRegion = x@snps$SubRegion, burden.function = WSS, 
                        maf.threshold = 0.5, get.effect.size = FALSE, 
                        alpha = 0.05, cores = 10)
}
\arguments{
  \item{x}{ A bed matrix, only needed if \code{burden="CAST"} or \code{burden="WSS"}}
  \item{NullObject}{ A list returned from \code{NullObject.parameters} }
  \item{genomic.region}{ A factor containing the genomic region of each SNP, \code{x@snps$genomic.region} by default, for example the CADD regions}
  \item{SubRegion}{ A vector containing subregions within each \code{genomic.region}, \code{x@snps$SubRegion} by default, for example genomic categories}
  \item{burden.function}{A function to compute the genetic score, \code{WSS} by default.}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score. Set at 0.5 by default}
  \item{get.effect.size}{ TRUE/FALSE: whether to return effect sizes of the tested \code{genomic.region} (OR for categorical phenotypes, betas for continuous phenotypes) }
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
  \item{cores}{ How many cores to use, set at 10 by default. Only needed if \code{NullObject$pheno.type = "categorical"} }
}
\details{
  This function will return results from the regression of the phenotype on the genetic score(s) for each genomic region. Within each genomic region, a subscore will be computed for each SubRegion and one test will be performed for each genomic.region.

  If only two groups of individuals are present, a classical logistic regression is performed.
  If more than two groups of individuals are present, a non-ordinal multinomial regression is performed,
  comparing each group of individuals to the reference group indicated by the argument \code{ref.level} in \code{NullObject.parameters}.
  The choice of the reference group won't affect the p-values, but only the Odds Ratios.
  In both types of regression, the p-value is estimated using the Likelihood Ratio test and the function \code{burden.mlogit}.  
}
\value{
  A dataframe with one row per genomic region and two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  If \code{get.effect.size=TRUE}, a list is returned with the previous dataframe in \code{$Asso} and with \code{effect}, a list containing matrices with three columns:
  \item{OR}{ The OR value(s) associated to the subscores in the regression. If there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR}
  \item{l.upper}{ The upper bound of the confidence interval of each OR}
}


\seealso{ \code{\link{NullObject.parameters}}, \code{\link{burden.subscores}},  \code{\link{CAST}}, \code{\link{WSS}}}
\examples{
#Import data in a bed matrix
#x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

#Add population
#x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
#x <- select.inds(x, superpop=="EUR")
#x@ped$pop <- droplevels(x@ped$pop)

#Group variants within CADD regions and genomic categories
#x <- set.CADDregions(x)

#Filter of rare variants: only non-monomorphic variants with
#a MAF lower than 2.5%
#and with a adjusted CADD score greater than the median
#x1 <- filter.adjustedCADD(x, filter = "whole", maf.threshold = 0.025)

#run null model, using the 1000Genome population as "outcome"
#x1.H0 <- NullObject.parameters(pheno = x1@ped$pop, ref.level = "CEU",
#                               RVAT = "burden", pheno.type = "categorical")

#run burden test WSS
#res.subscores <- burden.subscores(x1, NullObject = x1.H0, burden = WSS, 
#                                  get.effect.size=TRUE, cores = 1)
#res.subscores$Asso # p-values
#res.subscores$effect #OR values

}

