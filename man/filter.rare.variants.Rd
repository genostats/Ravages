\name{filter.rare.variants}
\alias{filter.rare.variants}
\title{
  Rare variants filtering
}
\description{
  Filter rare variants based on a MAF threshold and a given number of SNP per genomic region
}
\usage{
filter.rare.variants(x, ref.level, filter=c("whole", "controls", "any"), 
                     maf.threshold=0.01, min.nb.snps, group)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{ref.level}{ The level corresponding to the controls group, only needed if \code{filter=="controls"} }
  \item{filter}{ On which group the filter will be applied }
  \item{maf.threshold}{ The MAF threshold used to define a rare variant}
  \item{min.nb.snps}{ The minimum number of snps needed to keep a genomic region}
  \item{group}{ A vector or factor indicating the group of each individual, if missing, \code{x@ped$pheno} is taken}
}
\value{
  A bed.matrix with filtered variants
}
\details{
  If \code{filter="whole"}, only the variants having a MAF lower than the threshold in the entire sample are kept.

  If \code{filter="controls"}, only the variants having a MAF lower than the threshold in the controls group are kept.

  If \code{filter="any"}, only the variants having a MAF lower than the threshold in any of the groups are kept.
}
\examples{
#GRR values according to the SKAT formula
GRR.del <- GRR.matrix(genes.maf=Kryukov, select.gene="R1",
                      GRR = "SKAT", GRR.multiplicative.factor=2)
  
#Simulation of one group of 1,000 controls and two groups of 500 cases, 
#50% of deleterious variants, 5 genomic regions are simulated.
x <- random.bed.matrix(genes.maf=Kryukov, size = c(1000, 500, 500), 
                       baseline = c(0.001, 0.001), select.gene="R1",
                       GRR.matrix.del = GRR.del, prop.del = 0.5, 
                       same.variant=FALSE, fixed.variant.prop=TRUE,
                       genetic.model = "multiplicative", replicates = 5)
                             
#Filter of rare variants: only non-monomorphic variants with 
#a MAF lower than 1% and regions with at least 5 variants are kept.
x.filter <- filter.rare.variants(x, filter = "any", 
                                 maf.threshold = 0.01, min.nb.snps = 5)
}
