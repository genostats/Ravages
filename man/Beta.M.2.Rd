\name{Beta.M.2}
\alias{Beta.M.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Beta-M statistic
}
\description{
Calculates the Beta-M statistic and its p-value using a bed.matrix
}
\usage{
Beta.M.2(x, group, genomic_region, which.snps = rep(TRUE, ncol(x)), target = 0, B.max = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A bed.matrix
}
  \item{group}{
A factor containing the group of each individual
}
  \item{genomic_region}{
A vector containing the genomic region of each SNP
}
  \item{which.snps}{
A logical vector indicated which SNP to include. By default, all the SNP of the considered genomic region are included
}
  \item{target}{
The maximum number of times the observed statistic can be exceeded by a permutated statistic
}
  \item{B.max}{
The maximum number of permutations to perform
}
}
\details{
For each permutation, the group of each individual is sampled and a new statistic is calculated based on these data (permutated statistic). These values are then compared to the observed statistic calculated on the real data. To obtain the p-value, we divide the number of times when a permutated statistic is greater than the observed statistic (let name this value A) by the total number of permutations performed.
To limit the number of permutations, we can decide to stop the program when a certain value of A is reached given by the \code{target} argument. We limit as well the total number of permutations (\code{B.max}). To obtain the p-value, we divide \code{target} by the number of permutations needed to reach the \code{target} value. When the \code{target} value is not reached, we divide A by \code{B.max} and we use an approximation with the normal distribution to estimate the p-value.
}
\value{
A data frame containing the following fields:
  \item{genomic.region}{The name of the considered \code{genomic_region} }
  \item{stat}{The Beta-M statistic associated to each genomic region}
  \item{A}{The number of times the observed Beta-M statistic is exceeded by a permutated statistic}
  \item{B}{\code{B.max} if the \code{target} is not reached or the number of permutations used to reach the \code{target}}
  \item{p}{The p-value associated to the Beta-M test}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, group, genomic_region, which.snps = rep(TRUE, ncol(x)), 
    target = 0, B.max = 0) 
{
    r <- .Call("oz_beta_m", PACKAGE = "oz", x@bed, which.snps, 
        genomic_region, group, target, B.max)
    res <- data.frame(genomic.region = levels(genomic_region), 
        r)
    if (!is.null(res$A)) 
        res$p <- (res$A + 1)/(res$B + 1)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
