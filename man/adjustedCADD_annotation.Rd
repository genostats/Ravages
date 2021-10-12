\name{adjustedCADD.annotation}
\alias{adjustedCADD.annotation}
\title{
  Variant annotation with adjusted CADD scores
}
\description{
  Annotate variants with the adjusted CADD scores (CADD PHRED scores for coding, regulatory and intergenic regions)
}
\usage{
adjustedCADD.annotation(x, variant.scores, cores = 10, verbose = T)
}
\arguments{
  \item{x}{ A bed.matrix annotated with CADD regions using \code{set.CADDregions} }
  \item{variant.scores}{ A dataframe containing the adjusted CADD scores of the variants (Optional, useful to gain in computation time if the adjusted CADD scores of variants in the study are available)}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{verbose}{ Whether to display information about the function actions }
}
\value{
  The bed matrix x with adjusted CADD scores in \code{adjCADD}.
}
\details{
  Variants are directly annotated with the adjusted CADD scores in the function using the file "AdjustedCADD_v1.4_202108.tsv.gz" downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/ in the repository of the package Ravages or the scores of variants can be provided to \code{variant.scores} to gain in computation time (this file should contain 5 columns: the chromosome ('chr'), position ('pos'), reference allele ('A1'), alternative allele ('A2') and adjusted CADD scores ('adjCADD'). As CADD scores are only available for SNVs, only those ones will have an adjusted CADD score.  

  Those adjusted scores are used in the \code{RAVA.FIRST()} pipeline to filter rare variants.  
  
  As this function can take time when many variants are present, it is recommended to use this function chromosome by chromosome for large datasets or to fitler the bed matrix before the annotation.
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{RAVA.FIRST}}, \code{\link{filter.adjustedCADD}} 
}
\examples{
#Import 1000Genome data from region around LCT gene
#x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Annotate variants with adjusted CADD score
#x <- adjustedCADD.annotation(x)
}
