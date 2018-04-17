\name{group.mafs}
\alias{group.mafs}
\title{
  MAF computation for cases and controls
}
\description{
  Generates MAF values for cases and controls based on MAF in the general population and OR values
}
\usage{
  group.mafs(pop.maf, OR, baseline)
}
\arguments{
  \item{pop.maf}{ A vector of MAF in the general population}
  \item{OR}{ An OR matrix with one row per cases group and one column per variant }
  \item{baseline}{ A vector of baseline values, one for each cases group }
}
\value{
  A matrix of MAF values containing one column per variant and one row per group (the first one being the controls group)
}
\author{
Herve Perdry and Ozvan Bocher
}
\examples{
#Creation of an OR matrix with the same variants being deleterious, neutral and protective in the two groups of cases
OR <- OR.matrix.same.variant(length(Kryukov$maf[Kryukov$unit=="R1"], c(2,4), c(0.5,0.25), 0.2, 0.1)
#Creation of the MAF matrix
MAF <- group.mafs(pop.maf = Kryukov$maf[Kryukov$unit=="R1"], OR = OR, baseline = c(0.001,0.001))
}
