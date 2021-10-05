\name{set.CADDregions}
\alias{set.CADDregions}
\title{ Variants annotation based on 'CADD regions' and genomic categories }
\description{
Attributes CADD regions and genomic categories to variants based on their positions
}
\usage{
set.CADDregions(x, verbose = T)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{verbose}{ Whether to display information about the function actions }
}
\details{
  To attribute variants to CADD regions and genomic categories, the files "CADDRegions.2021.hg19.bed.gz" and "FunctionalAreas.hg19.bed.gz" will be downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/ in the repository of the package Ravages.  
  CADD regions are non-overlapping regions that have been defined in the whole genome to perform rare variant association tests in the \code{RAVA.FIRST()} pipeline.  
  It is recommended to use this function chromosome by chromosome for large datasets for time and memory managment.
}
\value{
The same bed matrix as x with three additional columns :
  \item{genomic.region}{ The CADD region of each variant}
  \item{SubRegion}{ The genomic category of each variant among 'Coding', 'Regulatory' or 'Intergenic'}
  \item{adjCADD.Median}{ The median of adjusted CADD of variants observed at least to times in GnomAD genomes r2.0.1}
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{RAVA.FIRST}}, \code{\link{filter.adjustedCADD}}, \code{\link{burden.subscores}}
}
\examples{
#Import 1000Genome data from region around LCT gene
#x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Group variants within CADD regions and genomic categories
#x <- set.CADDregions(x)
#table(x@snps$genomic.region) #CADD regions
#table(x@snps$SubRegion) #Genomic categories
}
