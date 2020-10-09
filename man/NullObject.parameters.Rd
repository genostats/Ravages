\name{NullObject.parameters}
\alias{NullObject.parameters}
\title{ Null Model for SKAT and burden tests }
\description{ Get the parameters under the null model to peforms burden tests or SKAT}
\usage{
  NullObject.parameters(pheno, RVAT, pheno.type, ref.level, data, formula)
}
\arguments{
  \item{pheno}{ The phenotype of each individual: a factor if \code{pheno.type = "categorial"}, and a numeric vector if \code{pheno.type = "continuous"} }
  \item{RVAT}{ The type of Rare Variant Association Test (RVAT) to perform: should be "burden" or "SKAT"}
  \item{pheno.type}{ The type of phenotype: "categorial" for binary or multinomial traits, or "continuous"}
  \item{ref.level}{ The reference group of individuals for the regression, only needed if \code{RVAT = "burden"} a,d \code{pheno.type = "categorial"}}
  \item{data}{ Optional, a matrix containing the covariates with one column per covariate and one row per individual }
  \item{formula}{ Optional, an R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
}
\details{
  This function gets the parameters under the null model for SKAT or the burden tests.  

  For burden tests, it computes the Log-Likelihood under the null model used to perform the Likelihood Ratio Test.   

  For SKAT, it computes the probabilites for each individual of belonging to each group based on the group sizes and the potential covariates.  

  If \code{formula} is missing, all columns from \code{data} will be included as covariates.
}
\value{
  A list containing different elements depending on the \code{RVAT} performed and the \code{pheno.type}.  
  See individual help pages for the corresponding functions:   

  - \code{\link{burden.NullObject}} if \code{RVAT = "burden"} and \code{pheno.type = "categorial"}   

  - \code{\link{burden.NullObject.continuous}} if \code{RVAT = "burden"} and \code{pheno.type = "continuous"}   

  - \code{\link{SKAT.NullObject}} if \code{RVAT = "SKAT"} and \code{pheno.type = "categorial"}  
 
  - \code{\link{SKAT.NullObject.continuous}} if \code{RVAT = "SKAT"} and \code{pheno.type = "continuous"}
}
\seealso{ \code{\link{burden.NullObject}}, \code{\link{burden.NullObject.continuous}}, \code{\link{SKAT.NullObject}}, \code{\link{SKAT.NullObject.continuous}} }


