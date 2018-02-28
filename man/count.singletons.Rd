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
Herve Perdry and Ozvan Bocher
}
\examples{
#List of OR parameters
my.pars <- list(OR.del = c(4, 8), prob.del = 0.2, prob.pro = 0.05)
#Simulation of genotypes with 1 replicate for 400 controls and two groups of 200 case with the scenario 2
x <- random.bed.matrix(pop.maf = Kryukov$maf[Kryukov$unit=="R1"], size = c(400, 200, 200), baseline = c(0.001, 0.001), replicates = 1, OR.pars = my.pars, scenario = 2)
S <- count.singletons(x)
}
