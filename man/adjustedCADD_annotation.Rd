\name{adjustedCADD.annotation}
\alias{adjustedCADD.annotation}
\title{
  SNVs and Indels annotation with adjusted CADD scores
}
\description{
  Annotate SNVs and Indels with the adjusted CADD scores (CADD PHRED scores for coding, regulatory and intergenic regions)
}
\usage{
adjustedCADD.annotation(x, SNVs.scores = NULL, indels.scores = NULL,
			cores = 10, verbose = T, path.data)
}
\arguments{
  \item{x}{ A bed.matrix annotated with CADD regions using \code{set.CADDregions} }
  \item{SNVs.scores}{ A dataframe containing the ADJUSTED CADD scores of the SNVs (Optional, useful to gain in computation time if the adjusted CADD scores of variants in the study are available)}
  \item{indels.scores}{ A dataframe containing the CADD PHREDv1.4 scores of the indels - Compulsory if indels are present in \code{x}}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{verbose}{ Whether to display information about the function actions }
  \item{path.data}{ The repository where data for RAVA-FIRST are or will be downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/}
}
\value{
  The bed matrix x with adjusted CADD scores in \code{adjCADD}.
}
\details{
  This function calls \code{adjustedCADD.annotation.SNVs} and \code{adjustedCADD.annotation.indels}. See the help of those two functions for more details.  
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{adjustedCADD.annotation.SNVs}, \link{adjustedCADD.annotation.indels}, \link{RAVA.FIRST}}, \code{\link{filter.adjustedCADD}} 
}
\examples{
#Import 1000Genome data from region around LCT gene
#x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Annotate variants with adjusted CADD score
#x <- adjustedCADD.annotation(x)
}
