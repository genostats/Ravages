\name{region.by.pos}
\alias{region.by.pos}
\title{
Gene grouping based on SNP positions
}
\description{
Performs a gene grouping based on SNP pair distance
}
\usage{
region.by.pos(x, nb.groups)
}
\arguments{
  \item{x}{
A bed.matrix
}
  \item{nb.groups}{
The maximum number of groups to make
}
}
\details{
The distance between two successive SNP is calculated. Then, the maximum allowed distance between two SNP in one genomic region is raised until the number of groups is smaller than the \code{nb.groups}.
Note that the final number of groups will be equal or smaller than \code{nb.groups}
}
\value{
A data frame containing:
  \item{chr}{The chromosome of the SNP}
  \item{id}{The identifiant of the SNP}
  \item{pos}{The position of the SNP}
  \item{gpe}{The computed group of the SNP being a number between 1 and the final number of groups}
}
\author{
Herve Perdry and Ozvan Bocher
}
\examples{
pos <- region.by.pos(x, nb.groups = 20000) #Separation of the SNP in 20000 groups at the most
}
