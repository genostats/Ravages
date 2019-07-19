\name{SKAT}
\alias{SKAT}
\title{ Multi group SKAT test }
\description{ Peforms the multi group SKAT test }
\usage{
SKAT(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, Pi, 
     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
     perm.target = 100, perm.max = 5e4, debug=FALSE)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{group}{ A factor defining groups to be compared }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{Pi}{Matrix of probabilities of belonging to each group }
  \item{weights}{ A vector with the weight of each variant. By default, the weight of each variant is inversely proportionnal to its MAF}
  \item{maf.threshold}{ The MAF above which variants are eliminated (default is to keep all variants)}
  \item{perm.target}{ The number of permutations to reach to estimate the p-value and stop the permutation procedure. If not reached, \code{perm.max} permutations will be used}
  \item{perm.max}{ The maximum number of permutations to perform to estimate the p-value, will be used if \code{perm.target} is not reached}
  \item{debug}{ Whether to print details about the permutations (moments estimators, skewness, kurtosis, ...) in addition to the p-values, FALSE by default }
}
\value{
A data frame containing for each genomic region the statistic values, informations about the permutations, and 
the corresponding p-values.
}
\examples{
#Import 1000Genome data from region around LCT gene
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Group variants within known genes
x <- set.genomic.region(x)

#Filter of rare variants: only non-monomorphic variants with
#a MAF lower than 2.5%
#keeping only genomic regions with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

#run SKAT using the 1000 genome population as "outcome"
#The population is available in the object LCT.pop
SKAT(x1, LCT.pop)

#Simulation of covariates
#with different probabilities in GBR/CEU/FIN and IBS/TSI
covar <- data.frame(sex=c(sample(0:1, sum(table(LCT.pop)[c("CEU", "GBR", "FIN")]),TRUE,c(0.2,0.8)), 
                          sample(0:1, sum(table(LCT.pop)[c("TSI", "IBS")]),TRUE,c(0.8,0.2))), 
                    u=runif(length(LCT.pop)))

#Compute probabilities of each individual belonging to the different groups
Pi.matrix.LCT <- Pi.matrix(group=LCT.pop, data=covar, formula= ~ sex, ref.level="CEU")

#SKAT taking into account the covariate
SKAT(x1, group=LCT.pop, Pi=Pi.matrix.LCT)
}
