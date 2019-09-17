\name{CAST}
\alias{CAST}
\title{ Cohort Allelic Sum Test}
\description{ Calculates the CAST genetic score }
\usage{
CAST(x, genomic.region = x@snps$genomic.region, maf.threshold = 0.01)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{maf.threshold}{ The MAF used for the definition of a rare variant, set as 0.01 by default}
}
\value{
A matrix containing the CAST genetic score with one row per individual and one column per \code{genomic.region}
}
\seealso{ \code{\link{WSS}}, \code{\link{burden.mlogit}} }
\examples{
#Importation of data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)
#Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

# Group variants within known genes
x <- set.genomic.region(x)

# Filter variants with maf (computed on whole sample) < 0.025
# keeping only genomic region with at least 10 SNPs
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

# run burden test CAST, using the 1000Genome population as "outcome"
score.CAST <- CAST(x1, maf.threshold=0.025)
}
