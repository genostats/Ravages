\name{set.genomic.region.subregion}
\alias{set.genomic.region.subregion}
\title{ Variants annotation based on regions and subregions positions }
\description{
Attributes regions and subregions to variants based on given positions
}
\usage{
set.genomic.region.subregion(x, regions, subregions, split = TRUE)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{regions}{
    A dataframe in bed format (start is 0-based and end is 1-based) containing the regions with the fields : \code{Chr} (the chromosome of the gene),
    \code{Start} (the start position of the gene, 0-based), \code{End} (the end position of the gene, 1-based), and
    \code{Name} (the name of the gene - a factor),
  }
  \item{subregions}{ A dataframe containing the subregions in the same format as \code{regions} }
  \item{split}{ Whether to split variants attributed to multiple regions by duplicating this variants, set at TRUE by default}
}
\details{ 
Warnings: \code{regions$Name} and \code{subregions$Name} should be factors containing UNIQUE names of the regions, ORDERED in the genome order. 

If \code{x@snps$chr} is not a vector of integers, it should be a factor with same levels as \code{regions$Chr}.

If a variant is attributed to multiple genomic regions, it will be duplicated in the bed matrix with one row per genomic region if \code{split = TRUE}.  

This function can be applied before using \code{burden.subscores} to perform a functionally-informed burden tests with sub-scores for each \code{SubRegion} within each \code{genomic.region}. 
}
\value{
The same bed matrix as x with two additional columns: \code{x@snps$genomic.region} containing the annotation of the \code{regions} and \code{x@snps$SubRegion} containing the annotation of the \code{subregions}.
}
\seealso{
 \code{\link{set.genomic.region}}, \code{\link{burden.subscores}}
}
\examples{
#Import 1000Genome data from region around LCT gene
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)

#Group variants within known genes and 
#Within coding and regulatory regions
x <- set.genomic.region.subregion(x, 
 regions = genes.b37, subregions = subregions.LCT)
}
