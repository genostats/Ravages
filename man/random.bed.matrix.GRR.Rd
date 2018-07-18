\name{random.bed.matrix}
\alias{random.bed.matrix}
\title{
  Simulation of a bed.matrix
}
\description{
  Generates a bed.matrix with genotypes for cases and controls based on simulation parameters
}
\usage{
  random.bed.matrix.(pop.maf, size, baseline, replicates, GRR.pars, GRR.pars.2=NULL, OR.function, model=c("general", "additive", "dominant", "recessive"))
}
\arguments{
  \item{pop.maf}{ A vector containing the MAF in the general population for each variant }
  \item{size}{ A vector containing the size of the controls group and of each cases group }
  \item{baseline}{ A vector containing the baseline for each group of cases}
  \item{replicates}{ The number of simulations to perform }
  \item{GRR.pars}{ A list of parameters containing \code{OR.del}=a matrix of GRR values for deleterious variants with one row per group of cases and one column per variant, \code{OR.prot}=a vector of GRR values for protective variants with one row per group of cases and one column per variant (optional), \code{prob.del}=the probability for each variant of being deleterious and \code{prob.pro}=the probability for each variant of being protective. All these GRR values correspond to the GRR associated to the heterozygous genotype compared to the homozygous genotype for the reference allele.} 
  \item{GRR.pars.2}{ A list of parameters for the GRR associated to the homozygous genotype for the alternate allele, only needed if \code{model}="general". The format is the same as for the argument \code{GRR} but only \code{OR.del} and \code{OR.pro} are needed.}
  \item{OR.function}{ The function called to generate the OR matrix}
  \item{model}{ The genetic model of the disease }
}
\value{
  A bed.matrix with as much columns (variants) as replicates*length(pop.maf).
  The field x@snps$genomic.region contains the replicate number and the field x@ped$pheno contrains the group of each individual, "0" besing the controls group.
}
\details{
  If OR.matrix.fix is called, the deleterious, neutral and protective variants are not the same between the different groups of cases.
  If OR.matrix.same.fix.variant is called, the same variants are deleterious, neutral or protective in the different groups of cases.
  The genetic model of the disease needs to be specified in this function.
  If \code{model}="general", there is no link between the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele. Therefore, the user has to give two matrices of GRR, one for the heterozygous, the other for the homozygous for the alternate allele.
  If \code{model}="additive", we assume that the the GRR associated to the homozygous for the alternate allele is the square of the GRR associated to the heterozygous: only one GRR matrix is needed.
  If \code{model}="dominant", we assume that the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele are equal: only one GRR matrix is needed.
  If \code{model}="recessive", we assume that the GRR associated to the heterozygous is equal to 1: the GRR given is the one associated to the homozygous for the alternate allele.
  For the argument \code{pop.maf}, MAF from the file Kryukov or GnomADgenes available with the package Ravages can be used as MAF in the general population.
}
\author{
  Herve Perdry and Ozvan Bocher
}
\examples{
##To fulfill
}
