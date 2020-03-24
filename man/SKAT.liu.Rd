\name{SKAT.liu}
\alias{SKAT.liu}
\title{ Multi group SKAT test using Liu et al. approximation }
\description{ Peforms SKAT on two or more groups of individuals using Liu et al. approximation }
\usage{
SKAT(x, NullObject, genomic.region = x@snps$genomic.region,
     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5,
     useskew = TRUE)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{NullObject}{ A list returned from \code{SKAT.NullObject} }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{weights}{ A vector with the weight of each variant. By default, the weight of each variant is inversely proportionnal to its MAF, as it was computed in the original SKAT method}
  \item{maf.threshold}{ The MAF above which variants are removed (default is to keep all variants)}
  \item{useskwew}{ Whether to use the skewness for the chi-square approximation (TRUE), or the kurtosis (FALSE) }
}
\details{
  The method from Liu et al. 2008 is used where p-values are estimated using a chi-square approximation from moment'$

  If \code{useskew = FALSE}, the kurtosis is used instead of skewness in the chi-square approximation. This is equivalent to "liu.mod" in SKAT package.  
 
  This function is used bu \code{SKAT} when the sample size is larger than 2000.
}
\value{
A data frame containing for each genomic region:
  \item{stat}{ The observed statistics}
  \item{p.value}{ The p-value of the test}
}
\seealso{ \code{\link{SKAT.NullObject}}, \code{\link{SKAT.liu}}, \code{\link{SKAT.bootstrap}}}
\references{ Wu et al. 2011, \emph{Rare-variant association testing for sequencing data with the sequence kernel association test}, American Journal of Human Genetics \bold{82-93} doi:10.1016/j.ajhg.2011.05.029;
             Liu et al. 2008, \emph{A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables}, Computational Statistics & Data Analysis, doi:10.1016/j.csda.2008.11.025
}



