\name{compute.GRR.matrix}
\alias{compute.GRR.matrix}
\title{
  Compute a GRR matrix for the simulations
}
\description{
  Computes a GRR matrix based on simulation parameters for variants and genetic effects
}
\usage{
  compute.GRR.matrix(file.pop.maf = Kryukov, n.case.groups = 2, GRR = c("constant", "SKAT", "variable"), GRR.value = NULL, GRR.function = NULL, GRR.multiplicative.factor = 2, select.gene = levels(file.pop.maf$gene)[[1]])
}
\arguments{
  \item{file.pop.maf}{ A file containing the MAF in the general population (column maf) for variants within a gene (column gene), only needed if the GRR depends on the MAF of a variant. By default, the MAF from the file Kryukov are used}
  \item{n.case.groups}{ The number of cases groups, i.e. the number of groups where variants will have a GRR greater than 1, set as 2 by default}
  \item{GRR}{ How to calculate the GRR }
  \item{GRR.value}{ GRR value if all the variants have the same GRR. In this situation, \code{file.pop.maf} is not needed}
  \item{GRR.function}{ A function indicating how to calculate the GRR depending on the MAF in the general population, only needed if \code{GRR}="variable"}
  \item{GRR.multiplicative.factor}{ A vector containing the multiplicative factor for the GRR for each group of cases compared to the first group of cases}
  \item{select.gene}{ The name of the gene(s) to select in the file \code{file.pop.maf} if it contains multiple genes, "R1"  by default for the file Kryukov}
}
\value{
  A matrix containing the GRR for each variant (in column) in each group of cases (in lines)
}
\details{
  The GRR can be computed in three ways using the argument \code{GRR}.
  If \code{GRR}="constant", the same GRR is applied to all the variants with the value given to the argument \code{GRR.value}.
  If \code{GRR}="SKAT", the GRR are calculating using the formula from the paper presenting the SKAT method and thus depend on MAF.
  If \code{GRR}="variable", the GRR are calculating using a function given by the user depending only on the MAF in the general population.
  The GRR can be the same between the different groups of cases (\code{multiplicative.factor}=1 or different.There are \code{n.case.groups}-1 values in this vector.
  If the GRR depend on the MAF in the general population, a file needs to be passed to the argument \code{file.pop.maf}. This file should contain at least a column "gene" and a column "maf". By default, the file Kryukov is used, but the file GnomADgenes which contains MAF in the NFE population can also be used.
}
\author{
  Ozvan Bocher and Herve Perdry
}
\examples{
  #GRR calculated on the MAF from the first unit of the file Kryukov using the formula from the SKAT paper with two groups of cases,
  #the second group having GRR values twice as high as the first one
  GRR.del <- compute.GRR.matrix(GRR = "SKAT", file.pop.maf = Kryukov, n.case.groups = 2,
                                GRR.multiplicative.factor=2, select.gene = levels(file.pop.maf$gene)[[1]])
  GRR.del[,1:5]
}
