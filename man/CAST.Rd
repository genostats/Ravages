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
\seealso{ \code{\link{WSS}} }
\examples{
#Data simulation with GRR values according to the SKAT formula
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, 
                      n.case.groups = 2, select.gene="R1",
                      GRR.multiplicative.factor=2)
                                
#Simulation of one group of 1,000 controls and two groups of 500 cases, 
#50% of deleterious, 5 genomic regions are simulated.
x <- random.bed.matrix(genes.maf=Kryukov, size = c(1000, 500, 500),
                       baseline = c(0.001, 0.001), same.variant=FALSE,
                       fixed.variant.prop=TRUE, select.gene="R1",
                       GRR.matrix.del = GRR.del, prop.del = 0.5, 
                       genetic.model = "multiplicative", replicates=5)
                             
#Keep only rare variant
x <- filter.rare.variants(x, filter = "any", maf.threshold = 0.01, 
                          min.nb.snps = 5)
                            
#Compute the CAST score
score <- CAST(x, genomic.region = x@snps$genomic.region, 
              maf.threshold = 0.01)
}
