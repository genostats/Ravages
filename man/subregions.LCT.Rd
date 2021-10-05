\name{subregions.LCT}
\docType{data}
\alias{subregions.LCT}

\title{Exemple of functional categories}

\description{
  Example of arbitrary functional categories (coding or regulatory) in the LCT locus (bed format, GRCH37).  
  "Coding" corresponds to coding parts of the exons and "Regulatory" corresponds to everything that falls outside these coding regions.  

  Data contain the \code{Chr}, the \code{Start} position, the \code{End} position and the \code{Name} of all functional regions in the LCT locus.
}
  
\format{
  The data contain one dataframe with four columns:
  \describe{
    \item{\code{Chr}}{The chromosome of the gene}
    \item{\code{Start}}{The start position of the functional region (0-based)}
    \item{\code{End}}{The end position of the functional region (1-based)}
    \item{\code{Name}}{The name of the gene}
  }
}  

\seealso{ \code{\link{set.genomic.region.subregion}}, \code{\link{burden.subscores}} }
