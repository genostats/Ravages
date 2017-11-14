\name{mafs}
\alias{mafs}
\title{
Compute mafs from OR value
}
\description{
Computes mafs for cases and controls based on OR values
}
\usage{
mafs(maf, OR, baseline = 0.01)
}
\arguments{
  \item{maf}{
A vector of mafs, one for each SNP
}
  \item{OR}{
An OR value for all the SNPs or a vector of OR values, one for each SNP
}
  \item{baseline}{
The baseline risk of the disease in the genral population
}
}
\details{
}
\value{
A list with the fields:
  \item{maf.controls}{The mafs calculated for the group of controls}
  \item{maf.cases}{The mafs calculated for the group of cases}
}
\examples{
data(Kryukov2009European1800)
maf_controls <- mafs(Kru[,2], 1)$maf.controls
maf_cases <- mafs(Kru[,2], sample(c(1,3), length(Kru[,2], replace=TRUE, prob=c(0.8,0.2))$maf.cases #Each variant as a probability of 0.2 to have a deleterious OR of 3

}
