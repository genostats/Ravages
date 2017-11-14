\name{pos_gpe}
\alias{pos_gpe}
\title{
Gene grouping based on SNP positions
}
\description{
Performs a gene grouping using a bed.matrix based on the SNP positions
}
\usage{
pos_gpe(x, nb_groups)
}
\arguments{
  \item{x}{
A bed.matrix
}
  \item{nb_groups}{
The maximum number of groups to make
}
}
\details{
The distance between two successive SNP is calculated. Then, the maximum allowed distance between two SNP in one genomic region is raised as long as the number of groups is greater than the \code{nb_groups}.
}
\value{
A data frame containing:
  \item{chr}{The chromosome of the SNP}
  \item{id}{The name of the SNP}
  \item{pos}{The position of the SNP}
  \item{gpe}{The computed group of the SNP being the a number between 1 and \code{nb_groups}}
}
\examples{
pos <- pos_gpe(x, 20000) #Separation of the SNP in 20000 groups at the most

}
