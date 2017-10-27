\name{mafs}
\alias{mafs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute mafs from OR value
}
\description{
Computes mafs for cases and controls based on OR values
}
\usage{
mafs(maf, OR, baseline = 0.01)
}
%- maybe also 'usage' for other objects documented here.
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
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the fields:
  \item{maf.controls}{The mafs calculated for the group of controls}
  \item{maf.cases}{The mafs calculated for the group of cases}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hervé Perdry
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
maf_controls <- mafs(Kru[,2], 1)$maf.controls
maf_cases <- mafs(Kru[,2], sample(c(1,3), length(Kru[,2], replace=TRUE, prob=c(0.8,0.2))$maf.cases #Each variant as a probability of 0.2 to have a deleterious OR of 3

## The function is currently defined as
function (maf, OR, baseline = 0.01) 
{
    if (length(OR) == 1 & length(maf) > 1) 
        OR <- rep(OR, length(maf))
    if (length(maf) == 1 & length(OR) > 1) 
        maf <- rep(maf, length(OR))
    alpha <- (OR - 1) * (1 - baseline)
    beta <- 1 + (OR - 1) * (baseline - maf)
    p.tem <- ifelse(OR == 1, maf, (sqrt(beta^2 + 4 * alpha * 
        maf) - beta)/(2 * alpha))
    uu <- OR * p.tem/(1 - p.tem)
    p.cas <- uu/(1 + uu)
    list(maf.controls = p.tem, maf.cases = p.cas)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
