\name{same.OR.matrix}
\alias{same.OR.matrix}
\title{
  Creation of OR matrix
}
\description{
  Generates an OR matrix based on OR values and proportions
}
\usage{
  same.OR.matrix(n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro)
}
\arguments{
  \item{n.variants}{ The number of variants }
  \item{OR.del}{ OR values for deleterious variants (a values or a vector of values)}
  \item{OR.pro}{ OR values for protective variants (as many OR.pro values as OR.del values)}
  \item{prob.del}{ The probability for a variant of being deleterious }
  \item{prob.pro}{ The probability for a variant of being protective }
}
\details{
  The OR matrix contains one row per (OR.del, OR.pro) couple values.
  The deleterious, neutral and protective variants are different in each row.
  Each variant has the same probability of being deleterious or protective.
}
\value{
  A matrix containing OR values for each variant with one variant per column and one OR couple values per row.
}
\author{
  Herve Perdry & Ozvan Bocher
}
\example{
same.OR.matrix(20 , c(2,2), c(0.5,0.5), 0.2, 0.1)
#Generates OR values for 20 variants
#The deleterious, neutral and protective variants are different between the two groups but have the same OR values
}

