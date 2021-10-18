\name{rbm.GRR.power}
\alias{rbm.GRR.power}
\title{
  Power of RVAT based on simulations and theoretical calculations (CAST) with GRR
}
\description{
  Computes the power of the tests CAST, WSS and SKAT based on simulations with GRR and based on theoretical calculations for CAST
}
\usage{
rbm.GRR.power(genes.maf = Kryukov, size = c(500, 500), prev = 0.01, 
              GRR.matrix.del, GRR.matrix.pro = NULL, 
              p.causal = 0.5, p.protect = 0, same.variant = FALSE, 
              genetic.model=c("multiplicative", "general", "dominant", "recessive"), 
              select.gene, alpha = 2.5e-6, selected.controls = TRUE, 
              power.type = c("simulations", "theoretical"), verbose = TRUE, 
              RVAT = c("CAST", "WSS", "SKAT"), 
              SKAT.method = c("permutations", "theoretical"),
              maf.threshold = 0.01, replicates = 1000, cores = 10)
}
\arguments{
  \item{genes.maf}{ A dataframe containing at least the MAF in the general population (column maf) for variants with their associated gene (column gene), by default the file \code{Kryukov} is used}
  \item{size}{ A vector containing the size of each group (the first one being the control group)}
  \item{prev}{ A vector containing the prevalence of each group of cases}
  \item{GRR.matrix.del}{ A list containing the GRR matrix associated to the heterozygous genotype compared to the homozygous reference genotype as if all variants are deleterious. An additional GRR matrix associated to the homozygous for the alternate allele is needed if \code{genetic.genetic.model="general"}}
  \item{GRR.matrix.pro}{ The same argument as \code{GRR.matrix.del} but for protective variants}
  \item{p.causal}{ The proportion of causal variants in cases}
  \item{p.protect}{ The proportion of protective variants in cases among causal variants}
  \item{same.variant}{ TRUE/FALSE: whether the causal variants are the same in the different groups of cases}
  \item{genetic.model}{ The genetic model of the disease }
  \item{select.gene}{ Which gene to choose from \code{genes.maf$gene} if multiple genes are present. If missing, only the first level is kept.}
  \item{alpha}{ The significance level to compute the power}
  \item{selected.controls}{ Whether controls are selected controls (by default) or controls from the general population}
  \item{power.type}{ Whether to compute the power based on 'simulations' (by default) or 'theoretical' calculations (only for CAST)}
  \item{verbose}{ Whether to print details about the running function}
  \item{RVAT}{ On which RVAT among 'CAST', 'WSS' and 'SKAT' to compute power (only needed if \code{power.type="simulations"} }
  \item{SKAT.method}{ Which method to use to compute SKAT ppower, i.e. permutations or theoretical moments (cf \code{SKAT} documentation)}
  \item{maf.threshold}{ The maf threshold to consider a causal variant  (set at 0.01 by default)}
  \item{replicates}{ On how many replicates the power should be computed}
  \item{cores}{ How many cores to use for moments computation, set at 10 by default }
}
\value{
  A single value giving the power of CAST if \code{power.type="theoretical"} or the power of \code{RVAT} if \code{power.type="simulations"}. 
}
\details{
  Simulations are performed in the same was as in \code{rbm.GRR}. Please refer to the documentation of this function.  

  Theoretical power is only available for CAST for which a non-central Chi-squared is used.   

}
\seealso{ \code{\link{GRR.matrix}}, \code{\link{Kryukov}}, \code{\link{GnomADgenes}}, \code{\link{rbm.GRR}} }
\examples{
#GRR values calculated with the SKAT formula
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, 
                      n.case.groups = 2, select.gene = "R1",
                      GRR.multiplicative.factor=2)
                              
#Simulation of one group of 1,000 controls and two groups of 500 cases, 
#each one with a prevalence of 0.001
#with 50% of causal variants, 5 genomic regions are simulated.
rbm.GRR.power(genes.maf = Kryukov, size = c(1000, 500, 500), 
              prev = c(0.001, 0.001), GRR.matrix.del = GRR.del, 
              p.causal = 0.5, p.protect = 0, select.gene="R1",
              same.variant = FALSE, genetic.model = "multiplicative", 
              power.type="theoretical", cores = 1, alpha = c(0.001,2.5e-6))
}
