\name{group.mafs.GRR}
\alias{group.mafs.GRR}
\title{
  Calculates the MAF in the controls group and in each group of cases
}
\description{
  Generates MAF values for the controls group and each group of cases based on MAF in the general population and GRR values
}
\usage{
  group.mafs.GRR(pop.maf, GRR, GRR.2 = NULL, baseline, model = c("general", "additive", "dominant", "recessive"))
}
\arguments{
  \item{pop.maf}{ A vector of MAF in the general population }
  \item{GRR}{ A matrix giving the GRR of the heterozygous compared to the homozygous for the reference allele with one row per cases group and one column per variant }
  \item{GRR.2}{ A matrix giving the GRR of the heterozygous compared to the homozygous for the reference allele with one row per cases group and one column per variant, only need if \code{model}="general"}
  \item{baseline}{ A vector containing the baseline values for each group of cases}
  \item{model}{ The genetic model of the disease }
}
\value{
  A matrix of MAF values with one column per variant and one row per group (the first one being the controls group)
}
\details{
  The genetic model of the disease needs to be specified in this function.
  If \code{model}="general", there is no link between the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele. Therefore, the user has to give two matrices of GRR, one for the heterozygous, the other for the homozygous for the alternate allele.
  If \code{model}="additive", we assume that the the GRR associated to the homozygous for the alternate allele is the square of the GRR associated to the heterozygous: only one GRR matrix is needed.
  If \code{model}="dominant", we assume that the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele are equal: only one GRR matrix is needed.
  If \code{model}="recessive", we assume that the GRR associated to the heterozygous is equal to 1: the GRR given is the one associated to the homozygous for the alternate allele.
  The \code{baseline} corresponds to the proportion of each sub-group of cases in the population. It is used only to calculate the MAF in the controls group.
  For the argument \code{pop.maf}, MAF from the file Kryukov or GnomADgenes available with the package Ravages can be used as MAF in the general population
}
\author{
  Herve Perdry and Ozvan Bocher
}
\examples{
##To fulfill
}
