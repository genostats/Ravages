\name{Kryukov2009European1800}
\docType{data}

\title{Kryukov data set}

\description{
  These data are the data from Kryukov et al, 2009 containing deep resequencing data for the European population.
  These data contain 200 units corresponding to 200 genes. For each unit, the data set contains the maf in the European population, the selection coefficient and the position of each variant.
}

\usage{ data(Kryukov2009European1800) }

\format{
  The data contains one dataframe called "Kryukov" with four columns:
  \describe{
    \item{unit}{The unit of each variant}
    \item{maf}{The maf of each variant in the European population}
    \item{selection.coefficient}{The selction coefficient of each variant in the European population}
    \item{position}{The position of each variant}
  }
}

\source{ The data were obtained from the SeqPower software (see also \url{http://www.bioinformatics.org/spower/input#data_download}). }

\references{Kryukov et al, 2009, \emph{Power of deep, all-exon resequencing for discovery of human trait genes}, Proceedings of the National Academy of Sciences, DOI:10.1073/pnas.0812824106}

