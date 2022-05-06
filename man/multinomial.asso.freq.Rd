\name{multinomial.asso.freq}
\alias{multinomial.asso.freq}
\title{
  Single variant association test with categorical phenotype
}
\description{
  Performs an association test between categorical phenotypes and single variants
}
\usage{
multinomial.asso.freq(x, pheno = x@ped$pheno, ref.level,
                      test = c("Genotypic", "Allelic"), get.effect.size = F,
                      min.maf.threshold = 0.05)
}
\arguments{
  \item{x}{ A bed matrix, only needed if \code{burden="CAST"} or \code{burden="WSS"}}
  \item{pheno}{ The phenotype of each individual: a factor if \code{pheno.type = "categorical"}, and a numeric vector if \code{pheno.type = "continuous"} }
  \item{ref.level}{ The reference group of individuals for the estimation of the effect size, only needed if \code{get.effect.size = T} }
  \item{test}{ Whether to perform the test on the three genotypes ("Genotypic") or on the two alleles ("Allelic")}
  \item{get.effect.size}{ TRUE/FALSE: whether to return effect sizes of the variants (OR) }
  \item{min.maf.threshold}{ MAF threshold used to define a frequent variant to apply single-variant test}
}
\value{
  A dataframe with one row per variant and three columns: the chromosome, position and p-value of each variant.  
  If \code{get.effect.size = T}, a list with \code{Asso} containing the previous dataframe and \code{OR} containing the OR in each group for each variant.
}
\details{
  This association test is based on a chi-square with the following number of df:  
  If \code{test = "Genotypic"}, (number of groups of individuals - 1)* 2  
  If \code{test = "Allelic"}, (number of groups of individuals - 1)  
}
\examples{
#Import data in a bed matrix
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

#Add population
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

#Perform association test
x.freq.asso <- multinomial.asso.freq(x, test = "Genotypic",
                                     pheno = x@ped$pop)
}
