\name{burden.weights}
\alias{burden.weights}
\title{
  Weights of variants for burden
}
\description{
  Computes weights of variants as performed in SKAT in the liability model
}
\usage{
  burden.weights(haplos, maf.threshold = 0.01)
}
\arguments{
  \item{haplos}{A matrix of haplotypes with one row per haplotype and one column per variant}
  \item{maf.threshold}{ The maf threshold to use for the definition of rare variants}
}
\value{
  A vector giving the weight of each variant
}
\details{
  The weight of each variant is computed using the same formula used in SKAT to build the liability model: w=-0.4*log10(MAF).
  Monomorphic variants and variants with a MAF lower than \code{maf.threshold} will have a weight of 0.
}
\examples{
#Load LCT matrix of haplotypes
data(LCT.haplotypes)
#Selection of LCT gene (136545910-136594712)
id.LCT <- subset(LCT.snps, pos>=136545910 & pos <= 136594712)$id
LCT.hap.LCTgene <- LCT.hap[,id.LCT]
weights.LCT <- burden.weights(LCT.hap.LCTgene)
}
