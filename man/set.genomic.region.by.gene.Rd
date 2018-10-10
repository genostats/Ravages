\name{set.genomic.region.by.gene}
\alias{set.genomic.region.by.gene}
\title{ Variants grouping based on gene positions }
\description{
Groups variants based on positions of known genes
}
\usage{
set.genomic.region.by.gene(x, genes, include.all = FALSE)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genes}{
    A dataframe containing the fields : \code{Chr} (the chromosome of the gene),
    \code{Start} (the start position of the gene), \code{End} (the end position of the gene), and
    \code{Gene_Name} (the name of the gene),
  }
  \item{include.all}{ TRUE/FALSE, whether or not to assigne all the variants into a genomic region. }
}
\details{ 
If \code{include.all == FALSE}, only the variants having their position between the \code{genes$Start} and the \code{genes$End} of a gene will be attributed to the corresponding gene.  
If \code{include.all == TRUE}, each variant will be included in the nearest gene.  
The files genes.b37 and genes.b38 can be used as \code{genes} argument
}
\value{
The same bed matrix as x with an additional column \code{x@snps$genomic.region} containing the group of each variant
}
\author{
Herve Perdry and Ozvan Bocher
}
