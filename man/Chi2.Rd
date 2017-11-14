\name{Chi2}
\alias{Chi2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Chi-squared test
}
\description{
Performs a chi-squared test and creates a vector using a matrix 
}
\usage{
Chi2(B_CAST, group)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B_CAST}{
A numeric matrix of 0 and 1 with one row per individual and one column per genomic region containing the CAST score
}
  \item{group}{
A factor containing the group of each individual
}
}
\details{
The classical chi-squared test has been modified to test if the expected value for each \code{genomic_region} is greater than 5 where a normal chi-squared test is performed. If the expected value is lower than 5, the p-values are computed by Monte Carlo simulation.
}
\value{
A vector containing for each \code{genomic_region} the p-value associated to the CAST burden test
}
\references{
%% ~put references to the literature/web site here ~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(gene_positions_hg37)
gpe_genes <- genes_gpe(x, hg37, include.all=TRUE)
B_CAST <- B_cast(x, gpe_genes$genes, maf.threshold=0.01)
pval <- apply(B, 2, Chi2, group=x@ped$centre) #Performs a chi-squared test between the different phenotypes in the bed matrix x for each genomic_region from the list gpe_genes 
}
