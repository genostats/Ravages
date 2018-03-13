\name{Beta.M}
\alias{Beta.M}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Beta-M statistic
}
\description{
Calculates the Beta-M statistic and its p-value
}
\usage{
Beta.M(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A bed.matrix}
  \item{group}{A factor containing the group of each individual}
  \item{genomic.region}{A vector containing the genomic region of each SNP}
  \item{which.snps}{A logical vector indicated which SNP to include. By default, all the SNP of the considered genomic region are included}
  \item{target}{The maximum number of times the observed statistic can be exceeded by a permutated statistic}
  \item{B.max}{The maximum number of permutations to perform}
}
\details{
For each permutation, the group of each individual is sampled and a new statistic is calculated based on these data (permutated statistic). These values are then compared to the observed statistic calculated on the real data. To obtain the p-value, we divide the number of times when a permutated statistic is greater than the observed statistic (let name this value A) by the total number of permutations performed.
To limit the number of permutations, we can decide to stop the program when a certain value of A is reached given by the \code{target} argument. We limit as well the total number of permutations (\code{B.max}). To obtain the p-value, we divide \code{target} by the number of permutations needed to reach the \code{target} value. When the \code{target} value is not reached, we divide A by \code{B.max} and we use an approximation with the normal distribution to estimate the p-value.
}
\value{
A data frame containing the following fields:
  \item{genomic.region}{The name of the considered \code{genomic.region} }
  \item{stat}{The Beta-M statistic associated to each genomic region}
  \item{nb.geq}{The number of times the observed statistic is equal or exceeded by a permutated statistic}
  \item{nb.eq}{The number of times the observed statistic is equal to a permutated statistic}
  \item{nb.perms}{\code{B.max} if the \code{target} is not reached or the number of permutations used to reach the \code{target}}
  \item{p.value}{The p-value calculated as (\code{nb.geq} - 0.5*\code{nb.eq} + 1) / (\code{nb.perms} + 1) }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Herve Perdry and Ozvan Bocher
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#List of OR parameters
my.pars <- list(OR.del = c(2, 4), prob.del = 0.2, prob.pro = 0.05)
#Simulation of genotypes with 10 replicates for 400 controls and two groups of 200 case with the same variants being deleterious, neutral or protective between the two groups of cases but different OR values.
x <- random.bed.matrix(pop.maf = Kryukov$maf[Kryukov$unit=="R1"], size = c(400, 200, 200), baseline = c(0.001, 0.001), replicates = 10, OR.pars = my.pars, OR.function = OR.matrix.same.variant)
#Keep only variants with MAF<1% in one of the three groups
x <- filter.rare.variants(x, filter = "any", maf.threshold = 0.01)
#Compute the Beta.M statistic on three groups
Beta.M(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e6)
#Compute the Beta.M statistic on two groups 
Beta.M(x, group = ifelse(x@ped$pheno==0, 0, 1), genomic.region = x@snps$genomic.region, which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e6)
}
