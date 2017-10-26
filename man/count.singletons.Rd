\name{count.singletons}
\alias{count.singletons}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count singletons
}
\description{
Calculates the number of singletons for each individual using a bed.matrix
}
\usage{
count.singletons(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A bed.matrix
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A vector giving the total number of singletons for each invidual
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
function (x) 
{
    w1 <- (x@snps$N1 == 1) & (x@snps$N2 == 0)
    w2 <- (x@snps$N1 == 1) & (x@snps$N0 == 0)
    C <- .Call("oz_count_alternative_alleles", PACKAGE = "oz", 
        x@bed, w1 | w2, w2)
    names(C) <- x@ped$id
    C
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
