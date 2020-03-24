\name{SKAT.NullObject}
\alias{SKAT.NullObject}
\title{ Null Model for multi group SKAT test }
\description{ Get the parameters under the null model to peforms SKAT}
\usage{
  SKAT.NullObject(group, data, formula)
}
\arguments{
  \item{group}{ A factor containing the group of each individual }
  \item{data}{ Optional, a matrix containing the covariates with one column per covariate and one row per individual }
  \item{formula}{ Optional, an R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
}
\details{
  This function computes the probabilites for each individual of belonging to each group based on the group sizes and the potential covariates.

  If \code{formula} is missing, all columns from data will be included as covariates.
}
\value{
  A list containing:
  \item{Pi.data}{ A matrix n.individuals x n.groups containing the probabilities that each individual belong to each group }
  \item{X}{ A matrix containing 1 in the first column for the intercept, and covariates from \code{data} and \code{formula}}
  \item{group}{ The factor containing the group of each individual}
}
\seealso{ \code{\link{SKAT}}, \code{\link{SKAT.liu}}, \code{\link{SKAT.bootstrap}} }

