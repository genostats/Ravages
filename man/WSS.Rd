\name{WSS}
\alias{WSS}
\title{ WSS genetic score}
\description{
  Caluclates the WSS genetic score
}
\usage{
WSS(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{group}{ A factor containing the group of each individual }
  \item{genomic.region}{ A vector containing the genomic region of each SNP }
}
\details{
}
\value{
A matrix containing the WSS genetic score with one row per individual and one column per \code{genomic.region}
}
\author{
Herve Perdry and Ozvan Bocher
}
\examples{
####Data simulation
#GRR values
GRR.del <- compute.GRR.matrix(GRR = "SKAT", file.pop.maf = Kryukov, n.case.groups = 2, GRR.multiplicative.factor=2, select.gene = "R1")
#Simulation of one group of 1,000 controls and two groups of 500 cases, each one with a prevalence of 0.001.
#50% of deleterious variants having a GRR calculated with the SKAT formula based on MAF in the file Kryukov are present and they are different between the two groups of cases.
#5 genomic regions are simulated.
x <- random.bed.matrix.GRR(file.pop.maf = Kryukov, size = c(1000, 500, 500), baseline = c(0.001, 0.001), GRR.matrix = GRR.del, prop.del = 0.5, prop.pro = 0, same.variant = FALSE, fixed.variant.prop = TRUE, genetic.model = "multiplicative", select.gene = "R1", replicates = 5)

#Compute the WSS score
WSS(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region )
}
