\name{genotypic.freq}
\alias{genotypic.freq}
\title{
  Genotypic frequencies per group of individuals
}
\description{
  Calculates the three genotypic frequencies in the controls group and each group of cases based on MAF in the general population and GRR values
}
\usage{
  genotypic.freq(file.pop.maf = Kryukov, GRR, GRR.2 = NULL, baseline, 
  				 genetic.model = c("general", "multiplicative", "dominant", "recessive"), select.gene=NULL)
}
\arguments{
  \item{file.pop.maf}{ A file containing the MAF in the general population (column maf) for variants grouped in genes (column gene), by default the file Kryukov from the package is used}
  \item{GRR}{ A matrix giving the GRR of the heterozygous genotype compared to the homozygous reference genotype with one row per cases group and one column per variant }
  \item{GRR.2}{ A matrix giving the GRR of the homozygous alternate genotype compared to the homozygous reference genotype with one row per cases group and one column per variant, only need if \code{genetic.model}="general"}
  \item{baseline}{ A vector containing the baseline values for each group of cases}
  \item{genetic.model}{ The genetic model of the disease }
  \item{select.gene}{ Which gene to choose in the column "gene" from \code{file.pop.maf} if multiple genes are present. If missing, the first level of \code{file.pop.maf$gene} is kept.}
}
\value{
  A matrix of MAF values with one column per variant and one row per group (the first one being the controls group)
}
\details{
  The genetic model of the disease needs to be specified in this function.
  If \code{genetic.model}="general", there is no link between the GRR associated to the heterozygous genotype and the GRR associated to the homozygous alternate genotype. Therefore, the user has to give two matrices of GRR, one for each of these genotypes.
  If \code{genetic.model}="multiplicative", we assume that the the GRR associated to the homozygous alternate genotype is the square of the GRR associated to the heterozygous genotype: only one GRR matrix is needed.
  If \code{genetic.model}="dominant", we assume that the GRR associated to the heterozygous genotype and the GRR associated to the homozygous alternate genotype are equal: only one GRR matrix is needed.
  If \code{genetic.model}="recessive", we assume that the GRR associated to the heterozygous genotype is equal to 1: the GRR given is the one associated to the homozygous alternate genotype.
  The \code{baseline} corresponds to the proportion of each sub-group of cases in the population. It is used only to calculate the MAF in the controls group.
  For the argument \code{file.pop.maf}, the file Kryukov or GnomADgenes available with the package Ravages can be used.
}
\author{
  Herve Perdry and Ozvan Bocher
}
\examples{
  #Construction of the GRR matrix using the formula from SKAT
  #GRR in the second group are twice as high as in the first group
  GRR.del <- compute.GRR.matrix(GRR = "SKAT", GRR.multiplicative.factor=2,
  								select.gene="R1")

  #Calculation of frequency in the three groups of individuals 
  #under a multilpicative model of the disease
  geno.freq.groups <- genotypic.freq(file.pop.maf = Kryukov, GRR = GRR.del, 
  									 baseline = c(0.001, 0.001), 
  									 genetic.model = "multiplicative", select.gene = "R1") 
}
