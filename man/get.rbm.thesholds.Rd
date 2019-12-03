\name{get.rbm.thresholds}
\alias{get.rbm.thresholds}
\title{
  Get simulations parameters for rbm.haplos.thresholds
}
\description{
  Get simulations parameters to use for genetic data simulations based on a libaility model   
  in \code{rbm.haplos.thresholds}.
}
\usage{
  get.rbm.thresholds(haplos, weights = -0.4*log10(colMeans(haplos)), maf.threshold, 
                     nb.causal, p.protect, h2, prev, normal.approx = TRUE)
}
\arguments{
  \item{haplos}{ A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{weights}{ A vector of weights for each variant to compute the burden used in the liability model. By default, wieghts=-0.4*log10(MAF)}
  \item{maf.threshold}{ The maf threshold to consider a rare variant  (set at 0.01 by default), variants with a MAF upper this threshold will have a weight of 0}
  \item{nb.causal}{ The number of causal variants}
  \item{p.protect}{ The proportion of protective variants among causal variants}
  \item{h2}{ The variance explained by the gene}
  \item{prev}{ A vector with the prevalence in each group of individuals}
  \item{normal.approx}{ TRUE/FALSE: whether to use the normal approximation to compute thresholds. Set at TRUE by default}
}
\value{
  A list with the following elements
  \item{burdens}{ A list containing the haplotypes burden in each group of individuals}
  \item{sd}{ A vector containing the standard deviation of non-genetic elements in each group of individuals}
  \item{thr1}{ A vector containing the lower threshold to use in each group of individuals}
  \item{thr2}{ A vector containing the upper threshold to use in each group of individuals}
}
\details{
  \code{nb.causal}, \code{p.protect}, \code{h2} and \code{prev} should be vectors of length corresponding to the number of groups to simulate. If they are of size 1, values will be duplicated; if they are longer than expected, vectors will be cut.  

  All monomorphic variants and variants with a MAF higher than \code{maf.threshold} will have a weight of 0. Causal variants are sampled among variants having weights greater than 0.  

  Burdens are then computed on causal variants weighted by \code{weights} for each haplotype in each group of individuals, and adjusted on the desired \code{h2}.  

  Finally, thresholds are chosen to respect the given \code{prev} in each group of individuals. These thresholds are computed using a normal approximation if \code{normal.approx=TRUE}, otherwise they are computed using a distribution from 1e6 sampled burdens. Be carreful if used with high \code{h2} values or low \code{prev} values.  

  To simulate a group of controls, \code{prev} needs to be set at 1, regardless of the other arguments.  
}
\seealso{ \code{\link{rbm.haplos.thresholds}} }
\examples{
  #Load LCT dataset for haplotype matrix
  data(LCT.haplotypes)

  #Parameters for simulations with one group of controls and two groups of cases
  #with the second group having half h2 and twice the prevalence
  #compared to the first one, 30 causal variants
  x.param <- get.rbm.thresholds(haplos=LCT.hap, maf.threshold = 0.01, nb.causal=30,
                                p.protect=0, h2=c(0.01, 0.01, 0.02), prev=c(1,0.01,0.005))
}

