\name{simulated.bedmatrix.haplo}
\alias{simulated.bedmatrix.haplo}
\title{
  Simulation of genetic data based on haplotypes
}
\description{
  Simulates genetic data with respect to allele frequency spectrum and linkage disequilibrium pattern observed on given haplotype data
}
\usage{
  simulated.bedmatrix.haplo(haplos, weights.variants=burden.weights, 
                            maf.threshold = 0.01, 
                            p.causal, p.protect = 0, nb.causal, 
                            h2, prev, normal.approx=TRUE, scenario, 
                            size=c(1000, 500, 500),
                            replicates=10, power, alpha.threshold=2.5e-6)
}
\arguments{
  \item{haplos}{ A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{weights.variants}{ A function to compute the weight of each variant}
  \item{maf.threshold}{ The maf threshold to filter rare variants for power calculation}
  \item{p.causal}{ The proportion of causal variants. Incompatible with \code{nb.causal}}
  \item{p.protect}{ The proportion of protective variants among causal variants}
  \item{nb.causal}{ The number of causal variants. Incompatible with \code{p.causal}}
  \item{h2}{ The heritability explained by the variants in the gene}
  \item{prev}{ A vector with the prevalence of each sub-phenotype, i.e with length=number of groups-1}
  \item{normal.approx}{ Whether to use a standard normal distribution to estimate the threshold for each sub-phenotype}
  \item{scenario}{ The function to use to simulate the data, i.e. the genetic scenario}
  \item{size}{ The sizes of each group of individuals}
  \item{replicates}{ The number of replicates to perform}
  \item{power}{ Needed only if power calulations are to be performed. A vector with at least "CAST", "WSS" or "SKAT" indicating for which test the power will be computed}
  \item{alpha.threshold}{ Needed only if power calulations are to be performed. The significance threshold to consider}
}
\value{
  A bed matrix with the simulated genetic data if no power calculation is performed
  A list with the bed matrix and a dataframe containing power levels for each test and each analysis  
}
\details{
  \code{prev} and \code{h2} should be vectors with length corresponding to the number of cases groups only in scenario "DVDR". 
  In all other scenarios, only one prevalence and one h2 value are needed.
}
\examples{
  #Subset of SKAT haplotypes matrix as example
  library(SKAT)
  data(SKAT.haplotypes)
  haplo.matrix <- SKAT.haplotypes$Haplotype[,1:20]
  #Prevalence of 1% in the two groups of cases, scenario DG
  x <- simulated.bedmatrix.haplo(haplos=haplo.matrix, nb.causal=2, prev=0.01, 
                                 h2=0.01, scenario="DG")

  #DVDR scenario with power calculation of the three tests 
  DVDR.withpower <- simulated.bedmatrix.haplo(haplos=haplo.matrix, nb.causal=2, 
                                              prev=c(0.01,0.005), h2=c(0.01, 0.02), scenario="DVDR",
                                            alpha.threshold=0.01, power=c("CAST", "WSS", "SKAT"))
  DVDR.withpower$x
  DVDR.withpower$power
}
