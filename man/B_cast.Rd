\name{B_cast}
\alias{B_cast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cast score
}
\description{
Calculates a CAST score for each individual and creates a matrix using a bed.matrix 
}
\usage{
B_cast(x, genomic_region, maf.threshold = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A bed.matrix
}
  \item{genomic_region}{
A vector containing the genomic region of each SNP
}
  \item{maf.threshold}{
The MAF used for the definition of a rare variant, set as 0.01 by default
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A numeric matrix of 0 and 1 with one row per individual and one column per \code{genomic_region} containing the CAST score
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hervé Perdry and Ozvan Bocher
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(gene_positions_hg37)
gpe_genes <- genes_gpe(x, hg37, include.all=TRUE)
B_CAST <- B_cast(x, gpe_genes$genes, maf.threshold=0.01)

## The function is currently defined as
function (x, genomic_region, maf.threshold = 0.01) 
{
    weights.0 <- ifelse((x@snps$maf == (1 - x@p) & x@snps$maf < 
        maf.threshold), 1, 0)
    weights.1 <- ifelse((x@snps$maf < maf.threshold), 1, 0)
    weights.2 <- ifelse((x@snps$maf == x@p & x@snps$maf < maf.threshold), 
        1, 0)
    if (length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | 
        length(weights.2) != ncol(x) | length(genomic_region) != 
        ncol(x)) {
        stop("x and weights dimensions mismatch")
    }
    if (!is.factor(genomic_region)) 
        stop("'genomic_region' is not a factor")
    B <- .Call("oz_burden2", PACKAGE = "oz", x@bed, nlevels(genomic_region), 
        genomic_region, weights.0, weights.1, weights.2)
    colnames(B) <- levels(genomic_region)
    rownames(B) <- x@ped$id
    Bcast <- (B > 0) + 0
    return(Bcast)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
