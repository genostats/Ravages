\name{Power}
\alias{Power}
\title{
  Power simulation of WSS, C.Alpha, Beta-M
}
\description{
  Calculates the power of WSS, C.Alpha, Beta-M based on simulated data
}
\usage{
  Power(alpha=0.05, filter=c("whole", "controls", "any"), maf.threshold=0.01, WSS=TRUE, C.ALPHA=TRUE, Beta.M=TRUE, model.pars)
}
\arguments{
  \item{alpha}{ Alpha threshold used to determine the significativity of the test}
  \item{filter}{ On which group the filter on rare variants will be applied }
  \item{maf.threshold}{ The MAF threshold used to define a rare variant}
  \item{WSS}{ Whether the WSS statistic is calculated}
  \item{C.ALPHA}{ Whether the C.ALPHA statistic is calculated}
  \item{Beta.M}{ Whether the Beta-M statistic is calculated}
  \item{model.pars}{A list of parameters with \code{pop.maf} = the MAF in the general population, \code{size} = a vector giving the size of each group of individuals, \code{baseline} = a vector giving the baseline in each group of cases, \code{replicates} = the number of replicates to calculate the power, \code{OR.pars} = a list of OR parameters, \code{scenario} = the scenario used to simulate the genotypes}
}
\value{
A data frame with the foloowing fields:
  \item{Power.WSS}{ The power of WSS in a non-pooled analysis }
  \item{Power.pooled.WSS}{ The power of WSS in a pooled analysis }
  \item{Power.C.ALPHA}{ The power of C.ALPHA in a non-pooled analysis }
  \item{Power.pooled.C.ALPHA}{ The power of C.ALPHA in a pooled analysis}
  \item{Power.Beta.M}{The power of Beta-M in a non-pooled analysis}
  \item{Power.pooled.Beta.M}{ The power of Beta-M in a pooled analysis}
  \item{Nb.replicates}{ The number of replicates (ie genomic regions) kept for each analysis}
}
\details{
  For each of the three statistics, two analysis are performed: a pooled and a non-pooled analysis.
  In the non-pooled analysis, each group of individuals is considered separately whereas in the pooled analysis, all the cases are considered as one group, which corresponds to a simple case-control association test.
  The power is calculated as the sum of tests having a p.value lower than \code{alpha} divided by the number of replicates kept in the analysis.
  This functions relies on the functions \code{random.bed.matrix.with.model} and \code{filter.rare.variants} used to generate the genotypes.
  OR.pars is a list of parameters containing \code{OR.del}=a vector of OR values for deleterious variants in each group of cases, \code{OR.prot}=a vector of OR values for protective variants in each group of cases, \code{prob.del}=the probability for each variant of being deleterious and \code{prob.pro}=the probability for each variant of being protective.
}
\author{
Herve Perdry and Ozvan Bocher
}
\examples{
#List of OR parameters with one group of cases having deleterious variants with OR=3 and the other group of cases OR=6
my.OR.pars <- list(OR.del=c(3,6), prob.del=0.2, prob.pro=0.05)
#List of model parameters with 400 controls and two groups of 200 cases, each one with a baseline of 0.001
my.pars <- list(Kryukov$maf[Kryukov$unit == "R1"], c(400,200,200), c(0.001,0.001), 1000, OR.pars=my.OR.pars, scenario=2)
#Power simulation on these parameters with the scenario 2 and the filter "any"
my.pow <- Power(alpha=0.05, filter="any", maf.threshold=0.01, WSS=TRUE, C.ALPHA=TRUE, Beta.M=TRUE, CAST=FALSE, SKAT=FALSE, model.pars=my.pars)
}
