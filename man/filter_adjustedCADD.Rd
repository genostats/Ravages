\name{filter.adjustedCADD}
\alias{filter.adjustedCADD}
\title{
  Variant filtering based on frequency and median adjusted CADD by CADD regions
}
\description{
  Filter rare variants based on a MAF threshold, a given number of SNP or a given cumulative MAF per genomic region and the median of adjusted CADD score for each CADD region
}
\usage{
filter.adjustedCADD(x, variant.scores = NULL, ref.level = NULL, 
                    filter=c("whole", "controls", "any"), 
                    maf.threshold=0.01, min.nb.snps = 2, 
                    min.cumulative.maf = NULL, 
                    group = NULL, cores = 10, verbose = T)
}
\arguments{
  \item{x}{ A bed.matrix annotated with CADD regions using \code{set.CADDregions} }
  \item{variant.scores}{ A dataframe containing the adjusted CADD scores of the variants (Optional, useful to gain in computation time if the adjusted CADD scores of variants in the study are available)}
  \item{ref.level}{ The level corresponding to the controls group, only needed if \code{filter=="controls"} }
  \item{filter}{ On which group the filter will be applied }
  \item{maf.threshold}{ The MAF threshold used to define a rare variant, set at 0.01 by default}
  \item{min.nb.snps}{ The minimum number of variants needed to keep a CADD region, set at 2 by default}
  \item{min.cumulative.maf}{ The minimum cumulative maf of variants needed to keep a CADD region}
  \item{group}{ A factor indicating the group of each individual, only needed if \code{filter = "controls"} or \code{filter = "any"}. If missing, \code{x@ped$pheno} is taken}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{verbose}{ Whether to display information about the function actions }
}
\value{
  A bed.matrix with filtered variants
}
\details{
  Variants are directly annotated with the adjusted CADD scores in the function using the file "AdjustedCADD_v1.4_202108.tsv.gz" downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/ in the repository of the package Ravages or the scores of variants can be provided to \code{variant.scores} to gain in computation time (this file should contain 5 columns: the chromosome ('chr'), position ('pos'), reference allele ('A1'), alternative allele ('A2') and adjusted CADD scores ('adjCADD'). As CADD scores are only available for SNVs, only those ones will be kept in the analysis.  

  If a column 'adjCADD' is already present in \code{x@snps}, no annotation will be performed and filtering will be directly on this column.  
  
  To use this function, a factor 'genomic.region' corresponding to the CADD regions and a vector 'adjCADD.Median' should be present in the slot \code{x@snps}. To obtain those two, use the function \code{set.CADDregions}.  
  
  Only variants with an adjusted CADD score upper than the median value are kept in the analysis. It is the filtering strategy applied in the \code{RAVA.FIRST()} pipeline.     

  If \code{filter="whole"}, only the variants having a MAF lower than the threshold in the entire sample are kept.

  If \code{filter="controls"}, only the variants having a MAF lower than the threshold in the controls group are kept.

  If \code{filter="any"}, only the variants having a MAF lower than the threshold in any of the groups are kept.  
  
  It is recommended to use this function chromosome by chromosome for large datasets.
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
 \code{\link{RAVA.FIRST}}, \code{\link{set.CADDregions}}, \code{\link{burden.subscores}}, \code{\link{filter.rare.variants} }
}
\examples{
#Import 1000Genome data from region around LCT gene
#x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Group variants within CADD regions and genomic categories
#x <- set.CADDregions(x)

#Annotate variants with adjusted CADD score
#and filter on frequency and median
#x.median <- filter.adjustedCADD(x, maf.threshold = 0.025, 
#                                min.nb.snps = 2)
}
