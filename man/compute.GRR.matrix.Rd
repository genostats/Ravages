\name{compute.GRR.matrix}
\alias{compute.GRR.matrix}
\title{
  Compute a GRR matrix for the simulations
}
\description{
  Computes a GRR matrix based on simulation parameters for variants and genetic effects
}
\usage{
  compute.GRR.matrix(file.pop.maf, n.case.groups, GRR, GRR.value, GRR.formula, same.GRR, GRR.multiplicative.factor, select.gene)
}
\arguments{
  \item{file.pop.maf}{ A file containing the MAF in the general population (column maf) for variants within a gene (column gene), only needed if the GRR depends on the MAF of a variant. By default, the MAF from the file Kryukov are used}
  \item{n.case.groups}{ The number of cases group, i.e. the number of groups where variants will have a GRR greater than 1}
  \item{GRR}{ How to calculate the GRR }
  \item{GRR.value}{ GRR value if all the variants have the same GRR, in this situation, file.pop.maf is not needed}
  \item{GRR.formula}{ A function indicating how to calculate the GRR depending on the MAF in the general population, only needed if \code{GRR}="variable"}
  \item{same.GRR}{ TRUE/FALSE: Whether the GRR values are the same between the different groups of cases}
  \item{GRR.multiplicative.factor}{ A vector containing the multiplicative factor for the GRR between each group of cases compared to the first group of cases}
  \item{select.gene}{ The name of the gene(s) to select in the file \code{file.pop.maf} if it contains multiple genes}
}
\value{
  A matrix containing the GRR for each variant (in column) in each group of cases (in lines)
}
\details{
  The GRR can be computed in three ways using the argument \code{GRR}.
  If \code{GRR}="constant", the same GRR is applied to all the variants with the value given to the argument \code{GRR.value}.
  If \code{GRR}="SKAT", the GRR are calculating using the formula from the paper presenting the SKAT method and thus depend on MAF.
  If \code{GRR}="variable", the GRR are calculating using a function given by the user depending only on the MAF in the general population.
  The GRR can be the same between the different groups of cases (\code{same.GRR}=TRUE) or different. In this situation, the user needs to define the multiplicative factor in each group of cases compared to the first group of cases: \code{n.case.groups}-1 values need to be specified.
  If the GRR depend on the MAF in the general population, a file need to be passed to the argument \code{file.pop.maf}. This file should contain at least a column "gene" and a column "maf". By default, the file Kryukov is used, but the file GnomADgenes which contains MAF in the NFE population can also be used.
}
\author{
  Ozvan Bocher and Herve Perdry
}
\examples{
##To fulfill
}
