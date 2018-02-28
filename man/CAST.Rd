\name{Cast}
\alias{Cast}
\title{ Cast statistic }
\description{ Calculates the CAST statistic and its p-value using a bed.matrix }
\usage{
Cast(x, group = x@ped$pheno, genomic.region =  genomic.region, maf.threshold = 0.01)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{group}{ A factor containing the group of each individual }
  \item{genomic.region}{ A factor defining the genomic region of each SNP }
  \item{maf.threshold}{
The MAF used for the definition of a rare variant, set as 0.01 by default
}
}
\details{
}
\value{
A data frame containing the following fields:
  \item{genomic.region}{The name of the considered \code{genomic.region} }
  \item{stat}{The WSS statistic for each genomic region }
  \item{p.value}{ The p-value calculated with a kruskall test}
}
\author{
Herve Perdry and Ozvan Bocher
}     
\examples{
my.pars <- list(OR.del = c(2, 4), prob.del = 0.2, prob.pro = 0.05)
#Simulation of genotypes with 10 replicates for 400 controls and two groups of 200 case with the scenario 2
x <- random.bed.matrix(pop.maf = Kryukov$maf[Kryukov$unit=="R1"], size = c(400, 200, 200), baseline = c(0.001, 0.001), replicates = 10, OR.pars = my.pars, scenario = 2)
#Keep only variants with MAF<1% in one of the three groups
x <- filter.rare.variants(x, filter = "any", maf.threshold = 0.01)
#Compute the CAST score
CAST(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region )
}
