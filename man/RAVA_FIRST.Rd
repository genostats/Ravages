\name{RAVA.FIRST}
\alias{RAVA.FIRST}
\title{
  RAVA-FIRST: RAre Variant Association using Functionally-InfoRmed STeps
}
\description{
  Analyse rare variants using the RAVA-FIRST approach based on CADD scores to group and filter rare variants
}
\usage{
RAVA.FIRST(x, SNVs.scores = NULL, indels.scores = NULL, ref.level, 
           filter=c("whole", "controls", "any"), 
           maf.threshold=0.01, min.nb.snps = 2, 
           min.cumulative.maf = NULL, group = NULL, 
           cores = 10, burden = TRUE, H0.burden, burden.parameters, 
           SKAT = TRUE, H0.SKAT, SKAT.parameters, verbose = TRUE, path.data)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{SNVs.scores}{ A dataframe with the columns 'chr', 'pos', 'A1', 'A2' and 'adjCADD' containing the ADJUSTED CADD scores of the SNVs (Optional, useful to gain in computation time if the adjusted CADD scores of variants in the study are available)}
  \item{indels.scores}{ A dataframe with the columns 'chr', 'pos', 'A1', 'A2' and 'PHRED_1.4' containing the CADD PHREDv1.4 scores of the indels - Compulsory if indels are present in \code{x}}
  \item{ref.level}{ The level corresponding to the controls group, only needed if \code{filter=="controls"} }
  \item{filter}{ On which group the MAF filter will be applied }
  \item{maf.threshold}{ The MAF threshold used to define a rare variant, set at 0.01 by default}
  \item{min.nb.snps}{ The minimum number of variants needed to keep a CADD region, set at 2 by default}
  \item{min.cumulative.maf}{ The minimum cumulative maf of variants needed to keep a CADD region}
  \item{group}{ A factor indicating the group of each individual, only needed if \code{filter = "controls"} or \code{filter = "any"}. If missing, \code{x@ped$pheno} is taken}
  \item{cores}{ How many cores to use, set at 10 by default}
  \item{burden}{ Whether to compute the burden test}
  \item{H0.burden}{A list returned from \code{NullObject.parameters} with \code{RVAT="burden"} }
  \item{burden.parameters}{ A list containing the parameters to use by \code{burden.subscores} for the burden analysis ('burden.function' and 'get.effect.size')}
  \item{SKAT}{ Whether to compute SKAT}
  \item{H0.SKAT}{A list returned from \code{NullObject.parameters} with \code{RVAT="SKAT"} }
  \item{SKAT.parameters}{ A list containing the parameters to use by \code{SKAT} ('get.moments', 'estimation.pvalue', 'params.sampling', 'debug')}
  \item{verbose}{ Whether to display information about the function actions }
  \item{path.data}{ The repository where data for RAVA-FIRST are or will be downloaded from https://lysine.univ-brest.fr/RAVA-FIRST/}
}
\value{
  A list containing the results for the burden analysis ('burden') and the results for the SKAT analysis ('SKAT'), along with information about CADD regions (positions, type of genomic categories overlapped by each region and median of adjusted CADD scores).
}
\details{
  Rare variants are analysed using the 'RAVA-FIRST' strategy composed of three steps:  
  - Rare variants are grouped in 'CADD regions' defined from the CADD scores of variants observed in GnomAD.  
  - Rare variant are selected within each CADD region based on an adjusted CADD score using a region-specific threshold corresponding to the median of scores observed in GnomAD in each region.  
  - Burden analysis is performed by integrating sub-scores for the coding, regulatory and intergenic categories within each CADD region. For SKAT analysis, a test for each CADD region is performed.  

  RAVA.FIRST() is based on the functions \code{set.CADDregions}, \code{filter.adjustedCADD}, \code{burden.subscores} and \code{SKAT}. Please refer to these functions for more information.
  Especially, refer to the functions \code{burden.subscores} and \code{SKAT} to get more information about what is need in \code{burden.parameters} and \code{SKAT.parameters}.   

  bedtools need to be installed on the system if indels are present or if scores are not provided for SNVs.  

  It is recommended to use this function chromosome by chromosome for large datasets.
}
\source{https://lysine.univ-brest.fr/RAVA-FIRST/}
\seealso{
  \code{set.CADDregions}, \code{filter.adjustedCADD}, \code{burden.subscores}, \code{SKAT}
}
\examples{
##This example takes a few minutes to run. To run it quickly, we advise the user to decrease the number of variants or to increase the number of cores used
#Import 1000Genome data from region around LCT gene
#x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)

#Add population
#x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]

#Select EUR superpopulation
#x <- select.inds(x, superpop=="EUR")
#x@ped$pop <- droplevels(x@ped$pop)

#Remove the multi-allelic variants
#x <- select.snps(x, !grepl(A2, pattern = ",") & !grepl(A1, pattern = ","))

#Keep only the rare variants
#x <- select.snps(x, maf <= 0.01 & maf > 0)

#Select the indels and make a dataframe to be used for annotation
#x.vcf <- subset(x@snps, nchar(A1)>1 | nchar(A2)>1)[,c("chr", "pos", "id", "A1", "A2")]

#Export the file for annotation 
#write.table(x.vcf, file = "indels.toannotate.vcf", col.names = F, row.names = F, quote = F, sep = "\t")

###############
#Annotate the file in on CADD website (https://cadd.gs.washington.edu/score)
#Please be careful here about which CADD version to use for annotation, here is v1.4
###############

#Download the results file, unzip it and load it in R
#x.indels <- read.table("indels.annotated.tsv")
#colnames(x.indels) <- c("chr", "pos", "A1", "A2", "CADD1.4", "PHRED_1.4")

#Set the NULL model for burden testing
#H0.burden <- NullObject.parameters(pheno = x@ped$pop, ref.level = "CEU", RVAT = "burden", pheno.type = "categorical")

#Run RAVA-FIRST
#x.RAVAFIRST = RAVA.FIRST(x, indels.scores = x.indels, H0.burden = H0.burden, SKAT = F, cores = 1, path.data = "RAVA-FIRST/Package", maf.threshold = 1)

}
