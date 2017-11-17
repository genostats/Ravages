\name{hg37}
\docType{data}

\title{hg37 gene positions}

\description{
  These data have been downloaded from Biomart on the Ensembl website with the hg37 version.
  The data contain the \code{Start} position, the \code{End} position, the \code{Gene_Name} and the \code{Chr} of all the genes in chromosomes 1 to 22 representing 54848 genes.
}
  
\usage{ data(hg37) }

\format{
  The data contain one dataframe with four columns:
  \describe{
    \item{Start}{The start postion of the gene}
    \item{End}{The end poisition of the gene}
    \item{Gene_Name}{The name of the gene}
    \item{Chr}{The chromosome of the gene}
  }
}  

\source{The data were obtained from the Ensembl website (see \url{http://grch37.ensembl.org/biomart/martview/}).}

