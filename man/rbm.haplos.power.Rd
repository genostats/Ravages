\name{rbm.haplos.power}
\alias{rbm.haplos.power}
\title{
  Power of RVAT based on simulations with haplotypes
}
\description{
  Computes the power of the tests CAST, WSS and SKAT based on simulations with haplotypes
}
\usage{
  rbm.haplos.power(haplos, freqs, weights = -0.4*log10(colMeans(haplos)), 
                   maf.threshold = 0.01, nb.causal, p.protect = 0, 
                   h2, prev, normal.approx = TRUE, size, verbose = TRUE,
                   alpha = 2.5e-6, RVAT = c("CAST", "WSS", "SKAT"),
                   simus.haplos = c("freqs", "liability"), 
                   replicates = 1000, rep.by.causal = 50, cores = 10)
}
\arguments{
  \item{haplos}{ A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{freqs}{ A matrix of haplotypes frequencies in each group of individuals, only needed if \code{simus.haplos = "freqs"} }
  \item{weights}{ A vector of weights for each variant to compute the burden used in the liability model, only needed if \code{simus.haplos = "liability"}. By default, wieghts=-0.4*log10(MAF). It can also be a single value that will be given to all variants}
  \item{maf.threshold}{ The maf threshold to consider a rare variant  (set at 0.01 by default). Only variants with a MAF upper than this threshold will be kept to compute RVAT power. If \code{simus.haplos="liability"}, variants with a MAF upper this threshold will have a weight of 0}
  \item{nb.causal}{ The number of causal variants, only needed if \code{simus.haplos = "liability"}}
  \item{p.protect}{ The proportion of protective variants among causal variants, only needed if \code{simus.haplos = "liability"}}
  \item{h2}{ The variance explained by the gene, only needed if \code{simus.haplos = "liability"}}
  \item{prev}{ A vector with the prevalence in each group of individuals, only needed if \code{simus.haplos = "liability"}}
  \item{normal.approx}{ TRUE/FALSE: whether to use the normal approximation to compute thresholds. Set at TRUE by default, only needed if \code{simus.haplos = "liability"}}
  \item{size}{ The sizes of each group of individuals}
  \item{verbose}{ Whether to display information about the function actions }
  \item{alpha}{ The significance level to compute the power}
  \item{RVAT}{ On which RVAT among 'CAST', 'WSS' and 'SKAT' to compute power (only needed if \code{power.type="simulations"} }
  \item{simus.haplos}{ Which method to simulate the data, if \code{simus.haplos="freqs"}, {rbm.haplos.freqs()} is used, otherwise \code{rbm.haplos.thresholds()} is used.}
  \item{replicates}{ The number of simulations to perform to estimate the power}
  \item{rep.by.causal}{ The number of time causal variants will be sampled}
  \item{cores}{ How many cores to use for moments computation, set at 10 by default }
}
\value{
  Power values of \code{RVAT}
}
\details{
  Simulations are perfromed accordingly to \code{rbm.haplos.thresholds()} or \code{rbm.haplos.freqs()}. Please refer to the corresponding manuals for more details on the simulation procedures.
}
\examples{
\donttest{
  #Simulations of 5 groups of individuals with haplotypes frequencies
  #from the 5 EUR populations

  #Load LCT dataset for haplotype matrix
  data(LCT.haplotypes)
  #Haplotypes for the variants in the LCT gene in the EUR population
  LCT.gene.hap <- LCT.hap[which(LCT.sample$super.population=="EUR"),
                          which(LCT.snps$pos>=136545410 & LCT.snps$pos<=136594750)]

  #Individuals from EUR
  LCT.sample.EUR <- subset(LCT.sample, super.population=="EUR")
  #Matrix of haplotypic frequencies
  LCT.freqs <- sapply(unique(LCT.sample.EUR$population), function(z)
                      ifelse(LCT.sample.EUR$population==z,
                             1/table(LCT.sample.EUR$population)[z], 0))

  #Simulation of genetic data for five groups of 50 individuals
  rbm.haplos.power(haplos=LCT.gene.hap, freqs=LCT.freqs, size=rep(50,5), 
                   replicates=5, rep.by.causal = 5,  RVAT = "CAST", cores = 1)
}
}