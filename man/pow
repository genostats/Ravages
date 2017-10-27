\name{pow}
\alias{pow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Power calculation based on simulated data
}
\description{
Calculates the power based on simulated data (simulated mafs and simulated genotypes) using mafs in the general population and OR values for each SNP
}
\usage{
pow(maf, nb.replicate, nb.groups, nb.controls, nb.case, prot = FALSE, filter.maf, nb.prot = 0.05, nb.del = 0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maf}{
A vector with mafs values in the general population
}
  \item{nb.replicate}{
The number of replicates used for the power calculation
}
  \item{nb.groups}{
The number of CASES groups to simulate
}
  \item{nb.controls}{
The number of controls 
}
  \item{nb.case}{
The number of cases in EACH group of cases
}
  \item{prot}{
A logicial indicating whether to include protective variant or not, set as FALSE by default
}
  \item{filter.maf}{
Which group to apply a maf filter of 0.01
}
  \item{nb.prot}{
The number of protective variants to include (a number between 0 and 1)
}
  \item{nb.del}{
The number of deleterious variants to include (a number between 0 and 1)
}
}
\details{
Mafs and genotypes are simulated for each replicate. 
If the filter "sample" is chosen, the filter is applied on the maf calculated on the entire dataset. If the filter "any" is chosen, all the SNPs having a maf lower than 0.01 in any of the groups are kept in the analysis.
For each replicate, two analyses are performed for each test: a "pooled" analyses (-P) considering all the cases as one group and a "no-pooled" analyses considering the different groups of cases separately.
To calculate the power, the number of times a test has its p-value under 5% is divided by the number of replicates.
}
\value{
A dataframe containing the following fiels:
  \item{test}{The statistical test performed (The P suffix corresponding to the "pooled" analysis)
}
  \item{nb_gps}{The number of groups considered for the analysis (2 for the "pooled" analysis)
}
  \item{OR}{The OR value of the deleterious variants
}
  \item{pow}{The calculated power associated to each test
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hervé Perdry and Ozvan Bocher
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(Kryukov2009European1800)
res <- pow2(Kru[,2], nb.replicate=1000, nb.groups=2, nb.controls=400, nb.case=200, prot=TRUE, filter.maf="any", nb.prot=0.05, nb.del=0.20)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
