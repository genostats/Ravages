\name{GnomADgenes}
\docType{data}

\title{GnomADgenes dataset}

\description{
  This dataframe contains MAF values in the NFE population from the GnomAD database.
}

\usage{data(GnomADgenes)}

\format{
  GnomADgenes is a dataframe with four columns:
  \describe{
    \item{chr}{ The chromosome of the variant}
    \item{pos}{ The position of the variant}
    \item{consequence}{ The functionnal consequence of the variant predicted by VEP}
    \item{gene}{ The gene of the variant predicted by VEP}
    \item{maf}{ The MAF of the variant in the NFE population}
  }
}

\source{The data were obtained from the GnomAD website (see \url{http://gnomad.broadinstitute.org/}) and the VEP website (see \url{https://www.ensembl.org/info/docs/tools/vep/}).}


