\name{SKAT}
\alias{SKAT}
\title{ Multi group SKAT test }
\description{ Peforms the multi group SKAT test }
\usage{
SKAT(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, Pi, 
     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, perm.target = 100, perm.max = 1e4)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{group}{ A factor defining groups to be compared }
  \item{genomic.region}{ A factor defining the genomic region of each variant }
  \item{Pi}{Matrix of probabilities of belonging to each group }
  \item{maf.threshold}{ The MAF above which variants are eliminated (default is to keep all variants)}
}
\value{
A data frame containing for each genomic region the statistic values, informations about the permutations, and 
the corresponding p-values.
}
