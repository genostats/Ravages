\name{SKAT}
\alias{SKAT}
\title{ Multi group SKAT test }
\description{ Peforms SKAT on two or more groups of individuals }
\usage{
SKAT(x, NullObject, genomic.region = x@snps$genomic.region, 
     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
     method = "size.based", params.bootstrap)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{NullObject}{ A list returned from \code{SKAT.NullObject} }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{weights}{ A vector with the weight of each variant. By default, the weight of each variant is inversely proportionnal to its MAF, as it was computed in the original SKAT method}
  \item{maf.threshold}{ The MAF above which variants are removed (default is to keep all variants)}
  \item{method}{ Which method to use among "size.based", "bootstrap", "liu", "liu.kurtosis". By default "size.based" that will choose the method depending on sample size (see \code{details}) }
  \item{params.bootstrap}{ A list containing the elements "perm.target", "perm.max", "debug" if \code{method = "boostrap"} }
}
\details{
  If \code{method} is not specified, bootstrap will be used if the sample size is lower than 2000, and the "liu" method otherwised.
  
  For small samples, p-values estimation is based on bootstrap sampling and a sequential procedure: permutated statistics are computed and each one is compared to the observed statistics.
  This method requires the arguments \code{perm.target}, \code{perm.max} and \code{debug} that should be given a list to \code{params.bootstrap}.
  If \code{params.bootstrap} is not specified, perm.target will be set at 100, perm.max at 1e4, and debug at FALSE.
  The boostrap progam stops when either \code{perm.target} or \code{perm.max} is reached.  
  P-values are then computed using a mixed procedure:   

  if \code{perm.target} is reached, the p-value is computed as : \code{perm.target} divided by the number of permutations used to reach \code{perm.target};   

  if \code{perm.max} is reached, the SKAT small sample procedure is used, and p-values are approximated using a chi-square distributions based on statistics' moments 1, 2 and 4 computed from the permutated values.     

  If \code{debug}==TRUE, more informations about the estimated statistics moments are given.  

  For large samples (>2000), the method from Liu et al. 2008 is used where p-values are estimated using a chi-square approximation from moment's calculation of the observed matrices used in the statistics. This is equivalent to "liu" in SKAT package.  
  
  If \code{method = "liu.kurtosis"}, the kurtosis is used instead of skewness in the chi-square approximation. This is equivalent to "liu.mod" in SKAT package. 
}
\value{
A data frame containing for each genomic region at least: 
  \item{stat}{ The observed statistics}
  \item{p.value}{ The p-value of the test. If \code{method = "bootstrap"}, \code{p.perm} if \code{perm.target} is reached, \code{p.chi2} if \code{perm.max} is reached. }

If \code{method = "bootstrap"}, additional fields are present:
  \item{p.perm}{ The p-value computed by permutations: number of times permutated is greater than observed statistics divided by the total number of permutations performed}
  \item{p.chi2}{ The p-value computed by the chi-square approximation using the SKAT small sample procedure}

Futhermore, if \code{debug}==TRUE, other informations are given about the moments estimation:
  \item{nb.gep}{ The number of times a permutated statistics is equal or greater than the observed statistics \code{stat}}
  \item{nb.eq}{ The number of times a permutated statistics is equal to the observed statistics \code{stat}}
  \item{nb.perms}{ The total number of simulations performed }
  \item{stat.mean}{ The mean of the permutated statistics}
  \item{M2, M3, M4}{ The statistics' moments 2, 3 and 4}
  \item{stat.var}{ The variance of the permutated statistics}
  \item{stat.skew}{ The skweness of the permutated statistics}
  \item{stat.kurt}{ The kurtosis of the permutated statistics}
}
\seealso{ \code{\link{SKAT.NullObject}}, \code{\link{SKAT.liu}}, \code{\link{SKAT.bootstrap}}}
\references{ Wu et al. 2011, \emph{Rare-variant association testing for sequencing data with the sequence kernel association test}, American Journal of Human Genetics \bold{82-93} doi:10.1016/j.ajhg.2011.05.029;
             Lee et al. 2012, \emph{Optimal Unified Approach for Rare-Variant Association Testing with Application to Small-Sample Case-Control Whole-Exome Sequencing Studies}, American Journal of Human Genetics, doi:10.1016/j.ajhg.2012.06.007;
             Liu et al. 2008, \emph{A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables}, Computational Statistics & Data Analysis, doi:10.1016/j.csda.2008.11.025
}

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
SKAT(x1, x1@ped$pop, perm.max = 500)
}
