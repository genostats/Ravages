\name{genes.positions}
\docType{data}
\alias{genes.positions}
\alias{genes.b37}
\alias{genes.b38}

\title{Genes positions}

\description{
  Positions of human genes. These data were downloaded from Biomart on the Ensembl website with the GRCh37 and GRCh38 version. Only the genes present in GnomAD were kept.

  The data contain the \code{Chr}, the \code{Start} position, the \code{End} position and the \code{Gene_Name} of all the genes in chromosomes 1 to 22 representing 19375 and 18278 genes
  respectively.
}
  
\usage{data(genes.b37)
data(genes.b38) }

\format{
  The data contain one dataframe with four columns:
  \describe{
    \item{\code{Chr}}{The chromosome of the gene}
    \item{\code{Start}}{The start position of the gene}
    \item{\code{End}}{The end position of the gene}
    \item{\code{Gene_Name}}{The name of the gene}
  }
}  

\source{The data were obtained from the Ensembl website (see \url{http://grch37.ensembl.org/biomart/martview/} and \url{http://www.ensembl.org/biomart/martview/}).}

