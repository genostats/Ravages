\name{WSS}
\alias{WSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WSS score
}
\description{
Calculates a WSS score for each individual and creates a matrix using a bed.matrix 
}
\usage{
WSS(x, genomic_region)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A bed.matrix
}
  \item{genomic_region}{
A vector containing the genomic region of each SNP
}
}
\details{
This function gives the individual WSS score by genomic region but does not give the p-value associated to this score. To calculate this p-value, a Kruskall-Wallis test must be performed as stated in the example.
}
\value{
A numeric matrix with one row per individual and one column per genomic region containing the WSS score
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
S <- WSS(x, gpe_genes$gpe)
pval <- apply(S, 2, function(b)kruskall.test(b, x@ped$pheno)$p.value  #Performs a kruskall-wallis test between the different phenotypes in the bed matrix x for each genomic_region from the list gpe_genes 

## The function is currently defined as
function (x, genomic_region) 
{
    Q <- (2 * x@snps$N2 + x@snps$N1 + 1)/(2 * (x@snps$N0 + x@snps$N1 + 
        x@snps$N2) + 2)
    W <- sqrt((nrow(x) - x@snps$NAs) * Q * (1 - Q))
    weights.0 <- ifelse(Q > 0.5, 2/W, 0)
    weights.1 <- 1/W
    weights.2 <- ifelse(Q < 0.5, 2/W, 0)
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
    return(B)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
