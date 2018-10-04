\name{region.by.gene}
\alias{region.by.gene}
\title{ Gene grouping based on gene positions }
\description{
Performs a gene grouping based on the positions of known genes
}
\usage{
region.by.genes(x, genes, include.all = FALSE)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genes}{
    A list containing the fields : \code{Gene_Name} (the name of the gene),a
    \code{Start} (the start position of the gene), \code{End} (the end position of the gene) 
    and \code{Chr} (the chromosome of the gene)
  }
  \item{include.all}{ Logical, whether or not to assigne all the SNP into a genomic region. }
}
\details{ 
If \code{include.all == FALSE}, only the SNPs having their position between the 
\code{genes$Start} and the \code{genes$End} of a gene will be attributed to the 
corresponding gene. If \code{include.all == TRUE}, all the SNPs will be included to the nearest gene.}
\value{
A data frame containing:
  \item{chr}{The chromosome of the SNP}
  \item{id}{The name of the SNP}
  \item{pos}{The position of the SNP}
  \item{gpe}{The computed group of the SNP being the \code{genes$Gene_Name}}
}
\author{
Herve Perdry and Ozvan Bocher
}
%\data
%The files genes.b37.RData and genes.b38.RData can be used as \code{genes} argument
\examples{
data(genes.b37)
gpe_genes <- region.by.gene(x, genes = genes.b37, include.all = TRUE) #Attribution of all the SNP in the genes from the list genes.b37
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
