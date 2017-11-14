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

x_controls <- geno.simu.controls(maf.controls=maf_controls, nb.controls=200)

## The function is currently defined as
function (maf.controls, nb.controls) 
{
    x.controls <- matrix(rbinom(nb.controls * length(maf.controls), 
        2, maf.controls), byrow = TRUE, nrow = nb.controls)
    x.controls <- as.bed.matrix(x.controls)
    x.controls@ped$pheno <- 0
    x.controls@ped$famid <- x.controls@ped$id <- paste("T", seq(1, 
        nb.controls), sep = "")
    return(x.controls)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
