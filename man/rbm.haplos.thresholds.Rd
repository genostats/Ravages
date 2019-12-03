\name{rbm.haplos.thresholds}
\alias{rbm.haplos.thresholds}
\title{
  Simulation of genetic data based on haplotypes and a libaility model
}
\description{
  Simulates genetic data with respect to allele frequency spectrum and linkage disequilibrium pattern observed on given haplotype data under a libaility model
}
\usage{
  rbm.haplos.thresholds(haplos, burdens, sd, thr1, size, replicates)
}
\arguments{
  \item{haplos}{ A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{burdens}{ A list (length=number of groups of individuals) of vector containing each haplotype burden}
  \item{sd}{A vector containing the standard deviation of the non-genetic composant of the phenotype in each group of individuals}
  \item{thr1}{ A vector containing the lower thresholds from the liability distribution to define each group of individuals}
  \item{size}{ The sizes of each group of individuals}
  \item{replicates}{ The number of simulations to perform}
}
\value{
  \item{x}{ A bed matrix with simulated genotypes}
}
\details{
  \code{sd}, \code{thr1}, \code{thr2}, \code{size} should be vectors of the same size.  
  \code{sd} corresponds to the standard deviation of the non-genetic composant. It is calculate as sqrt(1 - h²) where h² is the phenotypic variance explained by genetic effects.
  Haplotypic probabilities will be computed according to the fact that the burden liability is between thr1 and +Inf in each corresponding group of individuals. Two haplotypes are then sampled in each individual based on these probabilities.  
  The phenotype will be stored in \code{@ped$pheno}, and the simulation number is \code{@snps$genomic.region}.
}
\seealso{ \code{\link{get.rbm.thresholds}} }
\examples{
  #Load LCT dataset for haplotype matrix
  data(LCT.haplotypes)

  #Parameters for simulations with one group of controls and two groups of cases
  #with the second group having half h2 and twice the prevalence
  #compared to the first one, 30 causal variants
  x.param <- get.rbm.thresholds(haplos=LCT.hap, maf.threshold = 0.01, nb.causal=30, 
                                p.protect=0, h2=c(0.01, 0.01, 0.02), prev=c(1, 0.01, 0.005))

  #Simulation of 100 controls, and two groups of 50 cases using these parameters
  x <- rbm.haplos.thresholds(haplos=LCT.hap, burdens=x.param$burdens, sd=x.param$sd, 
                             thr1=x.param$thr1, size=c(100, 50, 50), replicates=5)
}
