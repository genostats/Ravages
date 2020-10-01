\name{SKAT.NullObject.continuous}
\alias{SKAT.NullObject.continuous}
\title{ Null Model for multi group SKAT test }
\description{ Get the parameters under the null model to peforms SKAT}
\usage{
  SKAT.NullObject.continuous(pheno, data = NULL, formula)
}
\arguments{
  \item{pheno}{ A numeric vector containing the phenotype value for each individual }
  \item{data}{ Optional, a matrix containing the covariates with one column per covariate and one row per individual }
  \item{formula}{ Optional, an R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
}
\details{
  This function computes the probabilites for each individual of belonging to each group based on the group sizes and the potential covariates.

  If \code{formula} is missing, all columns from data will be included as covariates.
}
\value{
  A list containing:
  \item{ymp}{ A matrix n.individuals x 1 containing the (y - pi_hat) values, i.e. the residuals from the regression of the phenotype on the potential covariates }
  \item{X}{ A matrix containing 1 in the first column for the intercept, and covariates from \code{data} and \code{formula}}
  \item{pheno}{ The phenotype of each individual as given}
  \item{P1}{ The variance matrix of \code{ymp} }
}
\seealso{\code{\link{NullObject.parameters}}, \code{\link{SKAT}}, \code{\link{SKAT.continuous}} }

