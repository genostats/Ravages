\name{filter.rare.variants}
\alias{filter.rare.variants}
\title{
  Filter rare variants
}
\description{
  Filter rare variants based on a MAF threshold and a given number of SNP per genomic region
}
\usage{
  filter.rare.variants(x, filter=c("whole", "controls", "any"), 
  					   maf.threshold=0.01, min.nb.snps)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{filter}{ On which group the filter will be applied }
  \item{maf.threshold}{ The MAF threshold used to define a rare variant}
  \item{min.nb.snps}{ The minimum number of snps needed to keep a genomic region}
}
\value{
  A bed.matrix with the SNPs kept by the filter
}
\details{
  If the filter "whole" is applied, only the SNPs having a MAF lower than the threshold in the entire sample are kept.
  If the filter "controls" is applied, only the SNPs having a MAF lower than the threshold in the controls group are kept.
  If the filter "any" is applied, only the SNPs having a MAF lower than the threshold in any of the groups are kept.
}
\author{
Herve Perdry and Ozvan Bocher
}
\examples{
#GRR values
  GRR.del <- compute.GRR.matrix(GRR = "SKAT", file.pop.maf = Kryukov, n.case.groups = 2, 
							    GRR.multiplicative.factor=2, select.gene = "R1")
  #Simulation of one group of 1,000 controls and two groups of 500 cases, 
  #each one with a prevalence of 0.001.
  #50% of deleterious variants having a GRR calculated with the SKAT formula based on MAF in the file 
  #Kryukov are present and they are different between the two groups of cases.
  #5 genomic regions are simulated.
  x <- random.bed.matrix.GRR(file.pop.maf = Kryukov, size = c(1000, 500, 500), baseline = c(0.001, 0.001), 
						     GRR.matrix = GRR.del, prop.del = 0.5, prop.pro = 0, same.variant = FALSE, 
						     genetic.model = "multiplicative", select.gene = "R1", replicates = 5)
  x
  #Filter of rare variants: only non-monomorphic variants with a MAF lower than 1% 
  #and regions with at least 5 variants are kept.
  x.filter <- filter.rare.variants(x, filter = "any", maf.threshold = 0.01, min.nb.snps = 5)
  x.filter
}
