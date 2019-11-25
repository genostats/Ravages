\name{random.bed.matrix.haplos}
\alias{random.bed.matrix.haplos}
\title{
  Simulation of genetic data based on haplotypes
}
\description{
  Simulates genetic data with respect to allele frequency spectrum and linkage disequilibrium pattern observed on given haplotype data
}
\usage{
  random.bed.matrix.haplos(haplos, weights=burden.weights, maf.threshold = 0.01, 
                           p.causal, p.protect = 0, nb.causal, 
                           h2, prev, size, replicates=1)
}
\arguments{
  \item{haplos}{ A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{weights}{ A function to compute the weight of each variant for the burden used in the liability}
  \item{maf.threshold}{ The maf threshold to consider a rare variant}
  \item{p.causal}{ The proportion of causal variants. Incompatible with \code{nb.causal}}
  \item{p.protect}{ The proportion of protective variants among causal variants}
  \item{nb.causal}{ The number of causal variants. Incompatible with \code{p.causal}}
  \item{h2}{ The variance explained by the gene}
  \item{prev}{ A vector with the prevalence in each group of individuals}
  \item{size}{ The sizes of each group of individuals}
  \item{replicates}{ The number of simulations to perform}
}
\value{
  \item{x}{ A bed matrix with simulated genotypes}
}
\details{
  \code{p.causal}/\code{nb.causal}, \code{p.protect}, \code{h2}, \code{prev} and \code{size} should be vectors of length corresponding to the number of groups to simulate. If they are of size 1, values will be duplicated.  

  If no function is given to \code{weights}, \code{burden.weights} is used to weight variants according to their MAF: w=-0.4*|log10(MAF)|. All monomorphic variants and variants with a MAF higher than \code{maf.threshold} will have a weight of 0.

  To generate a group of controls, it is sufficient to set \code{prev} at 1, regardless of the other arguments.  

  The phenotype of each individual will be stored in \code{@ped$pheno}, and the replicate number in \code{@snps$genomic.region}.
}
\examples{
  #Load LCT dataset for haplotype matrix
  #Selection of LCT gene (positions from GnomAD) and EUR superpop
  data(LCT.haplotypes)
  id.LCT <- subset(LCT.snps, pos>=136545410 & pos <= 136594750)$id
  haplo.matrix <- LCT.hap[which(LCT.sample$super.population=="EUR"),id.LCT]

  #Simulation of 200 controls, 100 individuals from two sub-phenotypes
  #with the second group having half h2 and twice the prevalence
  #compared to the first one
  #30 causal variants, different between the two groups of cases
  x <- random.bed.matrix.haplos(haplos=haplo.matrix, nb.causal=30, h2=c(0.01, 0.01, 0.02),
                                prev=c(1, 0.01, 0.005), size=c(200,100,100), replicates=5)
}
