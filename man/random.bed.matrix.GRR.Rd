\name{random.bed.matrix.GRR}
\alias{random.bed.matrix.GRR}
\title{
  Simulation of a bed.matrix
}
\description{
  Generates a bed.matrix with genotypes for cases and controls based on simulation parameters
}
\usage{
  random.bed.matrix.GRR(file.pop.maf, size, baseline, replicates, GRR.matrix, GRR.matrix.pro=NULL, 
  					 prop.del = 0.5, prop.pro = 0, same.variant = c(FALSE,TRUE), 
  					 fixed.variant.prop = c(TRUE,FALSE), select.gene,
  					 genetic.model=c("general", "multiplicative", "dominant", "recessive")
}
\arguments{
  \item{file.pop.maf}{ A file containing the MAF in the general population (column maf) for variants within a gene (column gene), by default the file Kryukov from the package is used}
  \item{size}{ A vector containing the size of the controls group and of each cases group }
  \item{baseline}{ A vector containing the baseline for each group of cases}
  \item{replicates}{ The number of simulations to perform }
  \item{GRR.matrix}{ A list containing the GRR matrix associated to the heterozygous compared to the homozygous for the reference allele for deleterious variants. An additional GRR matrix associated to the homozygous for the alternate allele is needed if \code{genetic.model}="general"}
  \item{GRR.matrix.pro}{ The same argument as \code{GRR.matrix} but for protective variants}
  \item{prop.del}{ The proportion of deleterious variants in the cases}
  \item{prop.pro}{ The proportion of protective variants in the cases}
  \item{same.variant}{ TRUE/FALSE: whether the causal variants are the same in the different groups of cases}
  \item{fixed.variant.prop}{ TRUE/FALSE: whether the arguments \code{prop.del} and \code{prop.pro} correspond to the final proportion of deleterious and protective variants (TRUE) or to the probability for each variant of being deleterious or protective (FALSE)}
  \item{genetic.model}{ The genetic model of the disease }
  \item{select.gene}{ Which gene to choose in the column "gene" from \code{file.pop.maf} if multiple genes are present. If missing, the first level of \code{file.pop.maf$gene} is kept.}
}
\value{
  A bed.matrix with as much columns (variants) as replicates*length(file.pop.maf$maf).
  The field x@snps$genomic.region contains the replicate number and the field x@ped$pheno contrains the group of each individual, "0" being the controls group.
}
\details{
  The genetic model of the disease needs to be specified in this function.
  If \code{model}="general", there is no link between the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele. Therefore, the user has to give two matrices of GRR, one for the heterozygous, the other for the homozygous for the alternate allele.
  If \code{model}="multiplicative", we assume that the the GRR associated to the homozygous for the alternate allele is the square of the GRR associated to the heterozygous: only one GRR matrix is needed.
  If \code{model}="dominant", we assume that the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele are equal: only one GRR matrix is needed.
  If \code{model}="recessive", we assume that the GRR associated to the heterozygous is equal to 1: the GRR given is the one associated to the homozygous for the alternate allele.
  For the argument \code{file.pop.maf}, MAF from the file Kryukov or GnomADgenes available with the package Ravages can be used as MAF in the general population. By default, only the first unit of the file Kyrukov is used (argument \code{select.gene}).
}
\author{
  Herve Perdry and Ozvan Bocher
}
\examples{
  #GRR values
  GRR.del <- compute.GRR.matrix(GRR = "SKAT", file.pop.maf = Kryukov, n.case.groups = 2, 
  								GRR.multiplicative.factor=2, select.gene = "R1")
  #Simulation of one group of 1,000 controls and two groups of 500 cases, 
  #each one with a prevalence of 0.001.
  #50% of deleterious variants having a GRR calculated with the SKAT formula based on MAF 
  #in the file Kryukov are present and they are different between the two groups of cases. 
  #5 genomic regions are simulated.
  x <- random.bed.matrix.GRR(file.pop.maf = Kryukov, size = c(1000, 500, 500), 
  							 baseline = c(0.001, 0.001), GRR.matrix = GRR.del, 
  							 prop.del = 0.5, prop.pro = 0, select.gene="R1",
  							 same.variant = FALSE, fixed.variant.prop = TRUE, 
  							 genetic.model = "multiplicative", replicates = 5)
}
