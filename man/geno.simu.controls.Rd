\name{geno.simu.controls}
\alias{geno.simu.controls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of genotypes in a group of controls
}
\description{
Simulates the genotype of individuals based on maf values for controls
}
\usage{
geno.simu.controls(maf.controls, nb.controls)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maf.controls}{
A vector with mafs values in a group of controls, ie in the general population
}
  \item{nb.controls}{
The number of controls for the simulation of genotypes
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A bed.matrix with the phenotype set as 1
}
\examples{
data(Kryukov2009European1800)
maf_controls <- mafs(Kru[,2], 1)$maf.controls
maf_cases <- mafs(Kru[,2], sample(c(1,3), length(Kru[,2], replace=TRUE, prob=c(0.8,0.2))$maf.cases

x_controls <- geno.simu.controls(maf.controls=maf_controls, nb.controls=200)

}
