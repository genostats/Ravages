\name{pos_gpe}
\alias{pos_gpe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene grouping based on SNP positions
}
\description{
Performs a gene grouping using a bed.matrix based on the SNP positions
}
\usage{
pos_gpe(x, nb_groups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A bed.matrix
}
  \item{nb_groups}{
The maximum number of groups to make
}
}
\details{
The distance between two successive SNP is calculated. Then, the maximum allowed distance between two SNP in one genomic region is raised as long as the number of groups is greater than the \code{nb_groups}.
}
\value{
A data frame containing:
  \item{chr}{The chromosome of the SNP}
  \item{id}{The name of the SNP}
  \item{pos}{The position of the SNP}
  \item{gpe}{The computed group of the SNP being the a number between 1 and \code{nb_groups}}
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
pos <- pos_gpe(x, 20000) #Separation of the SNP in 20000 groups at the most

## The function is currently defined as
function (x, nb_groups) 
{
    if (any(x@snps$pos == 0)) {
        stop("position equal to 0")
    }
    dist <- c(abs(diff(x@snps$pos)), 0)
    pos_d <- data.frame(chr = x@snps$chr, id = x@snps$id, pos = x@snps$pos, 
        d = dist)
    pos_d <- pos_d[which(pos_d$chr != 23 & pos_d$chr != 24 & 
        pos_d$chr != 26), ]
    seuil <- 1000
    while (length(which(pos_d$d > seuil)) > nb_groups) {
        seuil <- seuil + 1000
    }
    diff_pos <- diff(c(0, which(pos_d$d > seuil), length(pos_d$d)))
    gpe <- as.factor(rep.int(1:length(diff_pos), diff_pos))
    gpe_geno <- data.frame(pos_d, gpe = gpe)
    for (i in 1:21) {
        if (tail(gpe_geno$gpe[which(gpe_geno$chr == i)], 1) == 
            head(gpe_geno$gpe[which(gpe_geno$chr == i + 1)], 
                1)) {
            warnings("SNPs from 2 different chromosomes in the same group")
        }
    }
    return(gpe_geno)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
