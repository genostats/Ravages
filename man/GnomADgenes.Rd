\name{GnomADgenes}
\docType{data}
\alias{GnomADgenes}

\title{GnomADgenes dataset}

\description{
  This data frame contains MAF values in the Non-Finnish European (NFE) population from the GnomAD database.
}

\usage{data(GnomADgenes)}

\format{
  GnomADgenes is a dataframe with five columns:
  \describe{
    \item{chr}{ The chromosome of the variant}
    \item{pos}{ The position of the variant}
    \item{consequence}{ The functionnal consequence of the variant predicted by Variant Effect Predictor (VEP)}
    \item{gene}{ The gene of the variant predicted by VEP}
    \item{maf}{ The MAF of the variant in the NFE population}
  }
}

\source{The data were obtained from the GnomAD website (see \url{http://gnomad.broadinstitute.org/}) and the VEP website (see \url{https://www.ensembl.org/info/docs/tools/vep/}).}


