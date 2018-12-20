\name{set.genomic.region}
\alias{set.genomic.region}
\title{ Variants annotation based on gene positions }
\description{
Annotates variants based on positions of known genes
}
\usage{
set.genomic.region(x, genes = genes.b37, flank.width = 0L)
}
\arguments{
  \item{x}{ A bed.matrix }
  \item{genes}{
    A dataframe containing the fields : \code{Chr} (the chromosome of the gene),
    \code{Start} (the start position of the gene), \code{End} (the end position of the gene), and
    \code{Gene_Name} (the name of the gene),
  }
  \item{flank.width}{ An integer. Width of the flanking regions downstream and upstream the genes. }
}
\details{ 
If \code{x@snps$chr} is not a vector of integers, it should be a factor with same levels as \code{genes$Chr}. In practice
for human data we provide two data sets of genes on autosomes, \code{genes.b37} and \code{genes.b38}

If \code{flank.width} is null, only the variants having their position between the \code{genes$Start} and the \code{genes$End} of a gene will be attributed to the corresponding gene.  
When two genes overlap, variants in the overlapping zone will be attributed to the second one (in the order given by the position of their starting point on the genome).

If \code{flank.width} is a positive number, variants \code{flank.width} downstream or upstream a gene will be annotated as if in the gene. You can use \code{flank.width = Inf}
to have each variant included in the nearest gene.  
}
\value{
The same bed matrix as x with an additional column \code{x@snps$genomic.region} containing the annotation of each variant
}
\seealso{
 \code{\link{genes.b37}}, \code{\link{genes.b38}}
}
\examples{
x <- as.bed.matrix(LCT.gen, LCT.fam, LCT.bim)
x <- set.genomic.region(x)
table(x@snps$genomic.region)
}