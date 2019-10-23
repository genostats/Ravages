\name{WSS}
\alias{WSS}
\title{ WSS genetic score}
\description{
Caluclates the WSS genetic score
}
\usage{
WSS(x, genomic.region = x@snps$genomic.region)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genomic.region}{ A vector containing the genomic region of each variant }
}
\value{
A matrix containing the WSS genetic score with one row per individual and one column per \code{genomic.region}
}
\seealso{ \code{\link{CAST}} }
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

# Compute burden score WSS
score.WSS <- WSS(x1)
}
