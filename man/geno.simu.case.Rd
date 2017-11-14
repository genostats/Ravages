\name{geno.simu.case}
\alias{geno.simu.case}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of genotypes in a group of cases
}
\description{
Simulates the genotype of individuals based on maf values for cases
}
\usage{
geno.simu.case(maf.case, nb.case)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maf.case}{
A vector with mafs values in group of cases
}
  \item{nb.case}{
The number of cases for the simulation of genotypes
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A bed.matrix with the phenotype set as 0
}
\examples{
data(Kryukov2009European1800)
maf_controls <- mafs(Kru[,2], 1)$maf.controls
maf_cases <- mafs(Kru[,2], sample(c(1,3), length(Kru[,2], replace=TRUE, prob=c(0.8,0.2))$maf.cases

x_cases <- geno.simu.case(maf.case=maf.cases, 200)

}
