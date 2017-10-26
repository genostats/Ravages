\name{genes_gpe}
\alias{pos_gpe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gene grouping
}
\description{
Performs a gene grouping using a bed.matrix
}
\usage{
genes_gpe(x, genes, include.all = FALSE)
pos_gpe(x, nb_groups)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A bed.matrix
}
  \item{genes}{
A list containing the fields : “Gene_Name”(the name of the gene), “Start” (the start position of the gene), “End” (the end of the gene) and “Chr”(the chromosome number)
}
  \item{include.all}{
A logical indicated whether to include all SNPs in the grouping or not. Set as FALSE by default
  \item{nb_groups}{
The maximum number of groups to make
}
}
}
\details{
Gene grouping using SNP positions(pos_gpe): The distance between two successive SNP is calculated. Then, the maximum allowed distance between two SNP in one genomic region is raised as long as the number of groups is greater than the \code{nb_groups}. 
Gene grouping using known genes(genes_gpe): If \code{include.all} == FALSE, only the SNPs having their position between the \code{genes$Start} and the \code{genes$End} of a gene will be attributed to the corresponding gene. If \code{include.all == TRUE}, all the SNPs will be included. New gene positions are calculated by doing the average between the \code{genes$End} and the \code{genes$Start} of the next one. 

}
\value{
A data frame containing:
  \item{chr}{The chromosome of the SNP}
  \item{id}{The name of the SNP}
  \item{pos}{The position of the SNP}
  \item{gpe}{The computed group of the SNP being the \code{genes$Gene_Name} for genes_gpe or being a number between 1 and \code{nb_groups} for pos_gpe}
}
\data
The files gene_positions_hg37.RData and gene_positions_hg38.RData can be used as the \code{genes} argument
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

data(gene_positions_hg37)
gpe_genes <- genes_gpe(x, hg37, include.all=TRUE) #Attribution of all the SNP in the genes from the list hg37

## The function is currently defined as
function (x, genes, include.all = FALSE) 
{
    if (any(table(genes$Gene_Name) != 1) == TRUE) {
        a <- unique(as.character(genes$Gene_Name))
        genes$Gene_Name <- factor(as.character(genes$Gene_Name), 
            levels = a)
        w <- names(which(table(genes$Gene_Name) == 1))
        genes <- subset(genes, genes$Gene_Name \%in\% w)
        genes$Gene_Name <- droplevels(genes$Gene_Name)
    }
    if (include.all == TRUE) {
        genes_c <- genes
        for (i in 2:nrow(genes)) {
            if (genes$Start[i] > genes$End[i - 1]) {
                genes_c$Start[i] <- mean(c(genes$End[i - 1], 
                  genes$Start[i]))
                genes_c$End[i - 1] <- genes_c$Start[i]
            }
        }
        genes <- genes_c
    }
    gpe <- data.frame(chr = x@snps$chr, pos = x@snps$pos, id = x@snps$id, 
        gpe = factor(rep(NA, nrow(x@snps)), levels = levels(genes$Gene_Name)))
    for (i in 1:nrow(genes)) {
        w <- (x@snps$chr == genes$Chr[i] & x@snps$pos >= genes$Start[i] & 
            x@snps$pos <= genes$End[i])
        gpe$gpe[w] <- genes$Gene_Name[i]
    }
    gpe$gpe <- droplevels(gpe$gpe)
    return(gpe)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
