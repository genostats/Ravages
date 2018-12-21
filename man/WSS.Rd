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
#Data simulation with GRR values according to the SKAT formula
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, 
                      n.case.groups = 2, select.gene = "R1",
                      GRR.multiplicative.factor=2)
#Simulation of one group of 1,000 controls and two groups of 500 cases,
#50% of deleterious, 5 genomic regions are simulated.
x <- random.bed.matrix(genes.maf=Kryukov, size = c(1000, 500, 500),
                       baseline = c(0.001, 0.001),
                       GRR.matrix.del = GRR.del, prop.del = 0.5,
                       genetic.model = "multiplicative", replicates=5)
                          
#Keep only rare variant
x <- filter.rare.variants(x, filter = "any", maf.threshold = 0.01, 
                          min.nb.snps = 5)

#Compute the WSS score
score <- WSS(x, genomic.region = x@snps$genomic.region )
}
