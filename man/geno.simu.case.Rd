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
maf_controls <- mafs(Kru[,2], 1)$maf.controls
maf_cases <- mafs(Kru[,2], sample(c(1,3), length(Kru[,2], replace=TRUE, prob=c(0.8,0.2))$maf.cases

x_cases <- geno.simu.case(maf.case=maf.cases, 200)

## The function is currently defined as
function (maf.case, nb.case) 
{
    x.case <- matrix(rbinom(nb.case * length(maf.case), 2, maf.case), 
        byrow = TRUE, nrow = nb.case)
    x.case <- as.bed.matrix(x.case)
    x.case@ped$pheno <- 1
    x.case@ped$id <- x.case@ped$famid <- paste("C", seq(1, nb.case), 
        sep = "")
    return(x.case)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
