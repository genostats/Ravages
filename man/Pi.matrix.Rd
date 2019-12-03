\name{Pi.matrix}
\alias{Pi.matrix}
\title{Covariated-based individual probabilities}
\description{Compute probability of belonging to a group of individuals based on covariates}
\usage{
Pi.matrix(group, data, formula, ref.level)
}
\arguments{
  \item{group}{ A vector giving the group of each individual}
  \item{data}{ A matrix containing the covariates with one column per covariate and one row per individual}
  \item{formula}{ An R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
  \item{ref.level}{ The reference group of individuals for the regression}
}
\details{
The \code{ref.level} argument will not have any impact on the computed probabilities but is required to run the regression model.  
The results from this function can directly be used as the matrix of probabilities \code{Pi} in \code{SKAT}.
}
\value{
A matrix containing adjusted probabilities for each individual of belonging to each group with one row per individual and one column per group.
}
\seealso{ \code{\link{SKAT}} }
\examples{
#Simulation of covariates (one binary variable and one continuous variable)
#Binary variable: different probability depending on the origin
set.seed(1)
covar <- data.frame(sex=c(sample(0:1, sum(table(LCT.pop)[c("CEU", "GBR", "FIN")]), 
                                      TRUE, c(0.2,0.8)),
                          sample(0:1, sum(table(LCT.pop)[c("TSI", "IBS")]),TRUE,c(0.8,0.2))),
                    u=runif(length(LCT.pop)))

#Probabilities depending on the binary covariate
Pi.matrix.LCT.sex <- Pi.matrix(group=LCT.pop, data=covar, formula= ~ sex, ref.level="CEU")

#Probabilities based on the two covariates
Pi.matrix.LCT <- Pi.matrix(group=LCT.pop, data=covar, ref.level="CEU")
}
