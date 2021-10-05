\name{burden.subscores}
\alias{burden.subscores}
\title{
  Linear, logistic or multinomial regression on a multiple genetic scores within a genomic region
}
\description{
  Performs burden tests with subscores in the regression on categorial or continuous phenotypes
}
\usage{
burden.subscores(x, NullObject, genomic.region = x@snps$genomic.region, 
                 SubRegion = x@snps$SubRegion, burden.function = WSS, 
                 maf.threshold = 0.5, get.effect.size = FALSE, 
                 alpha = 0.05, cores = 10, verbose = TRUE)
}
\arguments{
  \item{x}{ A bed matrix}
  \item{NullObject}{ A list returned from \code{NullObject.parameters} }
  \item{genomic.region}{ A factor containing the genomic region of each SNP, \code{x@snps$genomic.region} by default, for example the CADD regions}
  \item{SubRegion}{ A vector containing subregions within each \code{genomic.region}, \code{x@snps$SubRegion} by default, for example genomic categories}
  \item{burden.function}{A function to compute the genetic score, \code{WSS} by default.}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score. Set at 0.5 by default}
  \item{get.effect.size}{ TRUE/FALSE: whether to return effect sizes of the tested \code{genomic.region} (OR for categorial phenotypes, betas for continuous phenotypes) }
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
  \item{cores}{ How many cores to use, set at 10 by default. Only needed if \code{NullObject$pheno.type = "categorial"} }
  \item{verbose}{ Whether to display information about the function actions }
}
\details{
  This function will return results from the regression of the phenotype on the genetic score(s) for each genomic region. Within each genomic region, a subscore will be computed for each SubRegion and one test will be performed for each genomic.region.  

  When used after \code{set.CADDregions}, it will perform a test by CADD region with one subscore by genomic category (coding, regulatory, intergenic) as in the \code{RAVA.FIRST()} strategy.  

  If only two groups of individuals are present, a classical logistic regression is performed.
  If more than two groups of individuals are present, a non-ordinal multinomial regression is performed,
  comparing each group of individuals to the reference group indicated by the argument \code{ref.level} in \code{NullObject.parameters}.
  The choice of the reference group won't affect the p-values, but only the Odds Ratios.
  In both types of regression, the p-value is estimated using the Likelihood Ratio test and the function \code{burden.mlogit}.  

  If the phenotype is continuous, a linear regression is performed using the function \code{burden.continuous}.  
  
  The type of phenotype is determined from \code{NullObject$pheno.type}.  
}
\value{
  A dataframe with one row per genomic region and two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  If \code{get.effect.size=TRUE}, a list is returned with the previous dataframe in \code{$Asso} and with \code{effect}, a list containing matrices with three columns:
  \item{OR/beta}{ The OR/beta value(s) associated to the subscores in the regression. For categorial phenotypes, if there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR/beta}
  \item{l.upper}{ The upper bound of the confidence interval of each OR/beta}
}


\seealso{\code{\link{RAVA.FIRST}}, \code{\link{NullObject.parameters}}, \code{\link{burden.continuous.subscores}}, \code{\link{burden.mlogit.subscores}},  \code{\link{CAST}}, \code{\link{WSS}}}
\examples{
#Import 1000Genome data from region around LCT gene
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

#Group variants within known genes and
#Within coding and regulatory regions
x <- set.genomic.region.subregion(x, regions = genes.b37, 
                                  subregions = subregions.LCT)

#Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

#Keep only variants with a MAF lower than 1%
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.01)

#run null model, using the 1000Genome population as "outcome"
x1.H0 <- NullObject.parameters(pheno = x1@ped$pop, ref.level = "CEU",
                               RVAT = "burden", pheno.type = "categorial")

#run functionally-informed burden test WSS in LCT
burden.subscores(select.snps(x1, genomic.region == "LCT"), 
                 NullObject = x1.H0, burden.function = WSS, 
                 get.effect.size=FALSE, cores = 1)


####Using the RAVA-FIRST approach with CDD regions
#Group variants within CADD regions and genomic categories
#x <- set.CADDregions(x)

#Filter of rare variants: only non-monomorphic variants with
#a MAF lower than 2.5%
#and with a adjusted CADD score greater than the median
#x1 <- filter.adjustedCADD(x, filter = "whole", maf.threshold = 0.025)

#run functionally-informed burden test WSS
#burden.subscores(x1, NullObject = x1.H0, burden.function = WSS, 
#                 get.effect.size=FALSE, cores = 1)

}

