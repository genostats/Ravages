\name{filter.rare.variants}
\alias{filter.rare.variants}
\title{
  Filter rare variants
}
\description{
  Filter rare variants based on a MAF threshold
}
\usage{
  filter.rare.variants(x, filter=c("whole", "controls", "any"), maf.threshold=0.01, min.nb.snps)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{filter}{ On which group the filter will be applied }
  \item{maf.threshold}{ The MAF threshold used to define a rare variant}
  \item{min.nb.snps}{ The minimum number of snps needed to keep a genomic region}
}
\value{
  A bed.matrix with the genomic regions containing enough rare variants
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
#List of OR parameters
my.pars <- list(OR.del=c(2, 4), prob.del=0.2, prob.pro=0.05)
#Simulation of genotypes with 1000 replicates for 400 controls and two groups of case with the scenario 2
x <- random.bed.matrix.with.model(Kryukov$maf[Kryukov$unit=="R1"], c(400, 200, 200), c(0.001, 0.001), 1000, OR.pars=my.pars, scenario=2)
#Filter of rare variants based on the controls group, only the genomic regions with 5 variants are kept
x.filter <- filter.rare.variants(x, filter="controls", maf.threshold=0.01, min.nb.snps=5)
}
