\name{SKAT.theoretical}
\alias{SKAT.theoretical}
\title{ Multi group SKAT test using Liu et al. approximation }
\description{ Peforms SKAT on two or more groups of individuals using Liu et al. approximation }
\usage{
SKAT.theoretical(x, NullObject, genomic.region = x@snps$genomic.region,
                 weights = (1 - x@snps$maf)**24, maf.threshold = 0.5,
                 estimation.pvalue = "kurtosis", cores = 10, debug = FALSE )
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{NullObject}{ A list returned from \code{SKAT.NullObject} }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{weights}{ A vector with the weight of each variant. By default, the weight of each variant is inversely proportionnal to its MAF, as it was computed in the original SKAT method}
  \item{maf.threshold}{ The MAF above which variants are removed (default is to keep all variants)}
  \item{estimation.pvalue}{ Whether to use the skewness ("skewness") or the kurtosis ("kurtosis") for the chi-square approximation }
  \item{cores}{ How many cores to use for moments computation, set at 10 by default }
  \item{debug}{ Whether to return the mean, standard deviation, skewness and kurtosis of the statistics. Set at FALSE by default }
}
\details{
  The method from Liu et al. 2008 is used where p-values are estimated using a chi-square approximation from moment'$

  If \code{estimation.pvalue = "kurtosis"}, the kurtosis is used instead of skewness in the chi-square approximation. This is equivalent to "liu.mod" in SKAT package.  
 
  This function is used by \code{SKAT} when the sample size is larger than 2000.
}
\value{
A data frame containing for each genomic region:
  \item{stat}{ The observed statistics}
  \item{p.value}{ The p-value of the test}

If \code{debug = TRUE}, the mean, standard deviation, skewness and kurtosis used to compute the p-value are returned
}
\seealso{ \code{\link{SKAT.NullObject}}, \code{\link{SKAT.theoretical}}, \code{\link{SKAT.bootstrap}}}
\references{ Wu et al. 2011, \emph{Rare-variant association testing for sequencing data with the sequence kernel association test}, American Journal of Human Genetics \bold{82-93} doi:10.1016/j.ajhg.2011.05.029;
             Liu et al. 2008, \emph{A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables}, Computational Statistics & Data Analysis, doi:10.1016/j.csda.2008.11.025
}



