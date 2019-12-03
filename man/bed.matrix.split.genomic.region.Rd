\name{bed.matrix.split.genomic.region}
\alias{bed.matrix.split.genomic.region}
\title{ Bed matrix for variants associated to multiple genomic regions}
\description{
Creates a new bed matrix with variants associated to multiple genomic regions being duplicated
}
\usage{
bed.matrix.split.genomic.region(x, changeID=TRUE, genomic.region=NULL, 
                                split.pattern=",")
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{changeID}{ TRUE/FALSE: whether to change the variants ID by including the gene name}
  \item{genomic.region}{ A vector containing the genomic region of each variant }
  \item{split.pattern}{ The character separating the genomic regions}
}
\value{
A bed matrix with variants assigned to multiple genomic regions being duplicated and the corresponding genomic regions separated
}
\details{
  If \code{changeID=TRUE}, variants will have new IDs being CHR:POS:A1:A2:genomic.region.
  The genomic region(s) associated to each varaint should be in \code{x@snps$genomic.region} or given as a vector to \code{genomic.region}. If both are present, \code{genomic.region} is used.
}


