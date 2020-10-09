\name{burden.NullObject.continuous}
\alias{burden.NullObject.continuous}
\title{ Null Model for burden tests on continuous phenotypes}
\description{ Get the parameters under the null model to peforms burden tests on continuous phenotypes}
\usage{
  burden.NullObject.continuous(pheno, data, formula)
}
\arguments{
  \item{pheno}{ A numeric vector containing the phenotype value for each individual }
  \item{data}{ Optional, a matrix containing the covariates with one column per covariate and one row per individual }
  \item{formula}{ Optional, an R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
}
\details{
  This function gets the parameters under the null model for the burden tests: it creates an object containing parameters to run \code{\link{burden}} on a continuous phenotype.

  If \code{formula} is missing, all columns from data will be included as covariates.
}
\value{
  A list containing:
  \item{pheno}{A numeric vector containing the phenotype value for each individual as given}
  \item{covar.toinclude}{ Which covariates to include in the regression, depends on the argument \code{formula} }
  \item{data}{ The \code{data} argument containing covariates, NULL if it was missing}
}
\seealso{ \code{\link{NullObject.parameters}} }

