\name{WSS}
\alias{WSS}
\title{ WSS score }
\description{
Calculates a WSS score for each individual and creates a matrix using a bed.matrix 
}
\usage{
WSS(x, genomic_region)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genomic_region}{ A vector containing the genomic region of each SNP }
}
\details{
This function gives the individual WSS score by genomic region but does not give the p-value associated to this score. To calculate this p-value, a Kruskall-Wallis test must be performed as stated in the example.
}
\value{
A numeric matrix with one row per individual and one column per genomic region containing the WSS score
}
\examples{
data(gene_positions_hg37)
gpe_genes <- genes_gpe(x, hg37, include.all=TRUE)
S <- WSS(x, gpe_genes$gpe)
pval <- apply(S, 2, function(b)kruskall.test(b, x@ped$pheno)$p.value  #Performs a kruskall-wallis test between the different phenotypes in the bed matrix x for each genomic_region from the list gpe_genes 

}
