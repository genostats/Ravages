\name{genes.b38}
\docType{data}

\title{Gene positions in version GRCh38}

\description{
  These data have been downloaded from Biomart on the Ensembl website with the GRCh38 version, only the genes present in GnomAD are kept.
  The data contain the \code{Chr}, the \code{Start} position, the \code{End} position and the \code{Gene_Name} of all the genes in chromosomes 1 to 22 representing 18278 genes.
}

\usage{ data(genes.b38) }

\format{
The data contain one dataframe with four columns:
  \describe{
    \item{Chr}{The chromosome of the gene}
    \item{Start}{The start postion of the gene}
    \item{End}{The end poisition of the gene}
    \item{Gene_Name}{The name of the gene}
  }
}
 
\source{The data were obtained from the Ensembl website (see \url{http://www.ensembl.org/biomart/martview/}).}
                      
  
