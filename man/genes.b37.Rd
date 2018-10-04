\name{genes.b37}
\docType{data}
\alias{genes.b37}

\title{genes.b37 gene positions}

\description{
  These data have been downloaded from Biomart on the Ensembl website with the genes.b37 version.
  The data contain the \code{Start} position, the \code{End} position, the \code{Gene_Name} 
  and the \code{Chr} of all the genes in chromosomes 1 to 22 representing 19375 genes.
}
  
\usage{ data(genes.b37) }

\format{
  The data contain one dataframe with four columns:
  \describe{
    \item{\code{Start}}{The start position of the gene}
    \item{\code{End}}{The end position of the gene}
    \item{\code{Gene_Name}}{The name of the gene}
    \item{\code{Chr}}{The chromosome of the gene}
  }
}  

\source{The data were obtained from the Ensembl website (see \url{http://grch37.ensembl.org/biomart/martview/}).}

