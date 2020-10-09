\name{burden.NullObject}
\alias{burden.NullObject}
\title{ Null Model for burden tests on multinomial phenotypes}
\description{ Get the parameters under the null model to peforms burden tests on categorial phenotypes}
\usage{
  burden.NullObject(group, ref.level, data, formula)
}
\arguments{
  \item{group}{ A factor containing the group of each individual }
  \item{ref.level}{ The reference group of individuals for the regression}
  \item{data}{ Optional, a matrix containing the covariates with one column per covariate and one row per individual }
  \item{formula}{ Optional, an R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
}
\details{
  This function get the parameters under the null model for the burden tests: it computes the Log-Likelihood under the null model used to performe the Likelihood Ratio Test.

  If \code{formula} is missing, all columns from data will be included as covariates.
}
\value{
  A list containing:
  \item{group}{A factor containing the group of each individual as given}
  \item{ref.level}{ The reference group of individuals for the regression as given}
  \item{H0.LogLik}{ The Log-Likelihood of the null model}
  \item{covar.toinclude}{ Which covariates to include in the regression, depends on the argument \code{formula} }
  \item{data}{ The \code{data} argument containing covariates, NULL if it was missing}
}
\seealso{ \code{\link{NullObject.parameters}} }

