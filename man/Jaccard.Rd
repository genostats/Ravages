\name{Jaccard}
\alias{Jaccard}
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

\examples{
}
