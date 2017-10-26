\name{Jaccard}
\alias{Jaccard}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Jaccard index
}
\description{
Calculates the Jaccard index for each pair of individuals using a bed.matrix
}
\usage{
Jaccard(x, maf.threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A bed.matrix
}
  \item{maf.threshold}{
The MAF used for the definition of a rare variant, set as 0.01 by default
}
}
\details{
The individuals carrying no rare variants will have a null Jaccard index with all the individuals including themselves.
}
\value{
A squared matrix giving the Jaccard index for each pair of individuals
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, maf.threshold = 0.01) 
{
    if (maf.threshold > 0.5) 
        stop("MAF threshold can't be this high")
    w <- (x@snps$maf <= maf.threshold)
    inverse <- (x@p >= 1 - maf.threshold)
    Ja <- .Call("oz_jaccard", PACKAGE = "oz", x@bed, w, inverse)
    J <- Ja$A/Ja$B
    J[Ja$B == 0] <- 0
    rownames(J) <- colnames(J) <- x@ped$id
    J
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
