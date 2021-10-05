\name{burden.continuous.subscores}
\alias{burden.continuous.subscores}
\title{
  Linear regression on a multiple genetic scores within a genomic region
}
\description{
  Performs burden tests with subscores in the regression on continuous phenotypes
}
\usage{
burden.continuous.subscores(x, NullObject, genomic.region = x@snps$genomic.region, 
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
  \item{get.effect.size}{ TRUE/FALSE: whether to return effect sizes of the tested \code{genomic.region} (OR for categorial phenotypes, betas for continuous phenotypes) }
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
  \item{cores}{ How many cores to use, set at 10 by default. Only needed if \code{NullObject$pheno.type = "categorial"} }
}
\details{
  This function will return results from the regression of the phenotype on the genetic score(s) for each genomic region. Within each genomic region, a subscore will be computed for each SubRegion and one test will be performed for each genomic.region.
}
\value{
  A dataframe with one row per genomic region and two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  If \code{get.effect.size=TRUE}, a list is returned with the previous dataframe in \code{$Asso} and with \code{effect}, a list containing matrices with three columns:
  \item{beta}{ The beta value(s) associated to the subscores in the regression}
  \item{l.lower}{ The lower bound of the confidence interval of each beta}
  \item{l.upper}{ The upper bound of the confidence interval of each beta}
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

#Simulation of a covariate + Sex as a covariate
#sex <- x1@ped$sex
#set.seed(1) ; u <- runif(nrow(x1))
#covar <- cbind(sex, u)

#Null model with the covariate sex and a continuous phenotype
#x1.H0.covar <- NullObject.parameters(pheno = x1@ped$pheno <- rnorm(nrow(x1)),
#                                     RVAT = "burden", pheno.type = "continuous",
#                                     data = covar, formula = ~ sex)

#WSS test
#res.subscores <-burden.continuous.subscores(x1, NullObject = x1.H0.covar, 
#                                            burden = WSS, get.effect.size=TRUE, cores = 1)
#res.subscores$Asso # p-values
#res.subscores$effect #beta values

}

