\name{SKAT}
\alias{SKAT}
\title{ Multi group SKAT test }
\description{ Peforms SKAT on two or more groups of individuals }
\usage{
SKAT(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
     perm.target = 100, perm.max = 5e4, debug=FALSE)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{group}{ A factor with the group of each individual }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{weights}{ A vector with the weight of each variant. By default, the weight of each variant is inversely proportionnal to its MAF, as it was computed in the original SKAT method}
  \item{maf.threshold}{ The MAF above which variants are removed (default is to keep all variants)}
  \item{perm.target}{ The number of times to exceed the observed statistics. If not reached, \code{perm.max} permutations will be used}
  \item{perm.max}{ The maximum number of permutations to perform to estimate the p-value, will be used if \code{perm.target} is not reached}
  \item{debug}{ Whether to print details about the permutations (moments estimators, skewness, kurtosis, ...) in addition to the p-values, FALSE by default }
}
\details{
  P-values estimation is based on permutations and a sequential procedure: permutated statistics are computed and each one is compared to the observed statistics, the progamms stops when either \code{perm.target} or \code{perm.max} is reached.  
  P-values are then computed using a mixed procedure:   

  if \code{perm.target} is reached, the p-value is computed as : \code{perm.target} divided by the number of permutations used to reach \code{perm.target};   

  if \code{perm.max} is reached, the SKAT small sample procedure is used, and p-values are approximated using a chi-square distributions based on statistics' moments 1, 2 and 4 computed from the permutated values.     

  If \code{debug}==TRUE, more informations about the estimated statistics moments are given.  
}
\value{
A data frame containing for each genomic region: 
  \item{stat}{ The observed statistics}
  \item{p.perm}{ The p-value computed by permutations: number of times permutated is greater than observed statistics divided by the total number of permutations performed}
  \item{p.chi2}{ The p-value computed by the chi-square approximation using the SKAT small sample procedure}
  \item{p.value}{ \code{p.perm} if \code{perm.target} is reached, \code{p.chi2} if \code{perm.max} is reached. }

If \code{debug}==TRUE, other informations are given about the moments estimation:
  \item{nb.gep}{ The number of times a permutated statistics is equal or greater than the observed statistics \code{stat}}
  \item{nb.eq}{ The number of times a permutated statistics is equal to the observed statistics \code{stat}}
  \item{nb.perms}{ The total number of simulations performed }
  \item{stat.mean}{ The mean of the permutated statistics}
  \item{M2, M3, M4}{ The statistics' moments 2, 3 and 4}
  \item{stat.var}{ The variance of the permutated statistics}
  \item{stat.skew}{ The skweness of the permutated statistics}
  \item{stat.kurt}{ The kurtosis of the permutated statistics}
}
\seealso{ \code{\link{Pi.matrix}}}
\references{ Wu et al, 2011, \emph{Rare-variant association testing for sequencing data with the sequence kernel association test}, American Journal of Human Genetics \bold{82-93} doi:10.1016/j.ajhg.2011.05.029}

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

#run SKAT using the 1000 genome EUR populations as "outcome"
#The population is available in the object LCT.pop
#The maximum number of permutations used is 500, please increase
#this number for a more accurate estimation of the p-values
SKAT(x1, LCT.pop, perm.max = 500)
}
