\name{adjustedCADD.annotation.indels}
\alias{adjustedCADD.annotation.indels}
\title{
  Indels annotation with adjusted CADD scores
}
\description{
  Annotate Indels with the adjusted CADD scores (CADD PHRED scores for coding, regulatory and intergenic regions)
}
\usage{
adjustedCADD.annotation.indels(x, variant.scores = NULL, 
			       cores = 10, verbose = T, path.data)
}
\arguments{
  \item{x}{ A bed.matrix annotated with CADD regions using \code{set.CADDregions} }
  \item{variant.scores}{ A dataframe containing the CADD PHREDv1.4 scores of the indels}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{verbose}{ Whether to display information about the function actions }
  \item{path.data}{ The repository where data for RAVA-FIRST are or will be downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/}
}
\value{
  The bed matrix x with adjusted CADD scores in \code{adjCADD}.
}
\details{
  Indels are directly annotated with the adjusted CADD scores in the function using the file "AdjustedCADD_v1.4_202204_indels.tsv.gz" downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/ in the repository of the package Ravages.  

  The adjusted CADD scores in "AdjustedCADD_v1.4_202204_indels.tsv.gz" have been computed using a set of 48M indels already annotated in the CADD website. If indels not present in this set are to be annotated, they will be given the same adjusted score as the indel with the nearest PHRED score v1.4 provided in \code{variant.scores} which should contain the chromosome ('chr'), position ('pos'), reference allele ('A1'), alternative allele ('A2') and PHRED CADD scores v1.4 ('PHRED_1.4').  

  Those adjusted scores are used in the \code{RAVA.FIRST()} pipeline to filter rare variants.  
  
  As this function can take time when a large number of SNVs are present, it is recommended to use this function chromosome by chromosome for large datasets or to fitler the bed matrix before the annotation.
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{adjustedCADD.annotation}, \link{adjustedCADD.annotation.SNVs}, \link{RAVA.FIRST}}, \code{\link{filter.adjustedCADD}} 
}
