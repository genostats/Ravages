\name{B_cast}
\alias{B_cast}
\title{ Cast score }
\description{ Calculates a CAST score for each individual and creates a matrix using a bed.matrix }
\usage{
B_cast(x, genomic_region, maf.threshold = 0.01)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genomic_region}{ A factor defining the genomic region of each SNP }
  \item{maf.threshold}{
The MAF used for the definition of a rare variant, set as 0.01 by default
}
}
\details{
}
\value{
A numeric matrix of 0 and 1 with one row per individual and one column per \code{genomic_region} containing the CAST score
}
\examples{
# data(gene_positions_hg37)
# gpe_genes <- genes_gpe(x, hg37, include.all=TRUE)
# Lactase data from gaston
# --> pas satisfaisant pour l'exemple, trop de variants rares
# --> est-ce que ça n'est pas possible d'utiliser la fonction 
#     de simulations pour construire un exemple ?
data(LCT)
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
regions <- factor( c(rep("R1",304), rep("R2", 303)) )
B_CAST <- B_cast(x, regions, maf.threshold = 0.05)

}
