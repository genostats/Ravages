\name{random.bed.matrix}
\alias{random.bed.matrix}
\title{
  Creation of a bed.matrix
}
\description{
  Generates a bed.matrix with genotypes for cases and controls based on OR values
}
\usage{
  random.bed.matrix.(pop.maf, size, baseline, replicates, OR.pars, scenario = c(1,2))
}
\arguments{
  \item{pop.maf}{ A vector containing the MAF in the general population for each variant }
  \item{size}{ A vector containing the size of each group(cases and controls) }
  \item{baseline}{ A vector containing the baseline for each group of cases}
  \item{replicates}{ The number of times genotypes are computed }
  \item{OR.pars}{ A list of parameters containing \code{OR.del}=a vector of OR values for deleterious variants in each group of cases, \code{OR.prot}=a vector of OR values for protective variants in each group of cases, \code{prob.del}=the probability for each variant of being deleterious and \code{prob.pro}=the probability for each variant of being protective} 
  \item{scenario}{ The scenario used to generate the OR matrix}
}
\value{
  A bed.matrix with as much columns as replicates*length(pop.maf).
  The field x@snps$genomic.region contains the replicate number.
}
\details{
  In the scenario 1, the same OR values are given for deleterious and protective variants but these variants are not the same between the differents group of cases.
  In the scenario 2, the same variants are deleterious, neutral and protective in the different groups of cases but the OR values are different.
  This function relies on the function \code{OR.matrix} or \code{same.OR.matrix} depending on the scenario and on \code{group.mafs}.
}
\author{
  Herve Perdry and Ozvan Bocher
}
\examples{
#List of OR parameters
my.pars <- list(OR.del = c(2, 4), prob.del = 0.2, prob.pro = 0.05)
#Simulation of genotypes with 10 replicates for 400 controls and two groups of case with the scenario 2
x <- random.bed.matrix(pop.maf = Kryukov$maf[Kryukov$unit=="R1"], size = c(400, 200, 200), baseline = c(0.001, 0.001), replicates = 10, OR.pars = my.pars, scenario = 2)
#With the scenario 1
my.pars <- list(OR.del = c(2, 2), prob.del = 0.2, prob.pro = 0.05)
x <- random.bed.matrix(pop.maf = Kryukov$maf[Kryukov$unit=="R1"], size = c(400, 200, 200), baseline = c(0.001, 0.001), replicates = 10, OR.pars = my.OR.pars, scenario = 1)
}
