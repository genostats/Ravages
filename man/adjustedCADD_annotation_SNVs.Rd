\name{adjustedCADD.annotation.SNVs}
\alias{adjustedCADD.annotation.SNVs}
\title{
  SNVs annotation with adjusted CADD scores
}
\description{
  Annotate SNVs with the adjusted CADD scores (CADD PHRED scores for coding, regulatory and intergenic regions)
}
\usage{
adjustedCADD.annotation.SNVs(x, variant.scores = NULL, 
			     cores = 10, verbose = T, path.data)
}
\arguments{
  \item{x}{ A bed.matrix annotated with CADD regions using \code{set.CADDregions} }
  \item{variant.scores}{ A dataframe containing the ADJUSTED CADD scores of the SNVs (Optional, useful to gain in computation time if the adjusted CADD scores of variants in the study are available)}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{verbose}{ Whether to display information about the function actions }
  \item{path.data}{ The repository where data for RAVA-FIRST are or will be downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/}
}
\value{
  The bed matrix x with adjusted CADD scores in \code{adjCADD}.
}
\details{
  SNVs are directly annotated with the adjusted CADD scores in the function using the file "AdjustedCADD_v1.4_202108.tsv.gz" downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/ in the repository of the package Ravages or the scores of variants can be provided to \code{variant.scores} to gain in computation time (this file should contain 5 columns: the chromosome ('chr'), position ('pos'), reference allele ('A1'), alternative allele ('A2') and adjusted CADD scores ('adjCADD').  

  Those adjusted scores are used in the \code{RAVA.FIRST()} pipeline to filter rare variants.  
  
  As this function can take time when a large number of SNVs are present, it is recommended to use this function chromosome by chromosome for large datasets or to fitler the bed matrix before the annotation.
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{adjustedCADD.annotation}, \link{adjustedCADD.annotation.indels}, \link{RAVA.FIRST}}, \code{\link{filter.adjustedCADD}} 
}
\examples{
#Import 1000Genome data from region around LCT gene
#x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Annotate variants with adjusted CADD score
#x <- adjustedCADD.annotation.SNVs(x)
}
