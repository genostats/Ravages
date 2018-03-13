\name{OR.matrix.same.variant}
\alias{OR.matrix.same.variant}
\title{
  Creation of OR matrix
}
\description{
  Generates an OR matrix based on OR values and probabilities of deleterious and protective variants
}
\usage{
  OR.matrix.same.variant(n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro)
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
  The deleterious, neutral and protective variants are the same in each row, corresponding to each group of individuals.
  The probabilities of being deleterious or protective are the same for all variants.
}
\value{
  A matrix containing OR values for each variant with one variant per column and one OR couple values per row.
} 
\author{
  Herve Perdry and Ozvan Bocher
}
\examples{
#Generates OR values for 20 variants, the deleterious, neutral and protective variants are the same between the two groups but with different OR values.
OR.matrix.same.variant(n.variants = 20 , OR.del = c(2,4), OR.pro = c(0.5,0.25), prob.del = 0.2, prob.pro = 0.1)
}
