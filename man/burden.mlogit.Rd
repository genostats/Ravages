\name{burden.mlogit}
\alias{burden.mlogit}
\title{
  Logistic or multinomial regression on a genetic score
}
\description{
  Performs a logistical or a non-ordinal multinomial regression on a genetic score
}
\usage{
burden.mlogit(x, group = x@ped$pheno, 
              genomic.region = x@snps$genomic.region, 
              burden, maf.threshold = 0.01, ref.level, 
              formula=NULL, data=NULL,
              get.OR.value = FALSE, alpha = 0.05)
}
\arguments{
  \item{x}{ A bed matrix, only needed if \code{burden="CAST"} or \code{burden="WSS"}}
  \item{group}{ A factor containing the group of each individual, \code{x@ped$pheno} by default}
  \item{genomic.region}{ A factor containg the genomic region of each SNP, \code{x@snps$genomic.region} by default, only needed if \code{burden}="CAST" or \code{burden}="WSS"}
  \item{burden}{"CAST" or "WSS" to directly compute the CAST or the WSS genetic score; or a matrix with one row per individual and one column per \code{genomic.region} containing another genetic score.}
  \item{maf.threshold}{ The MAF threshold to use for the definition of a rare variant in the CAST score}
  \item{ref.level}{ The reference group of individuals for the regression}
  \item{formula}{ An R formula corresponding to the regression model indicating which covariates from \code{data} to include in the model if only some of them are to be included}
  \item{data}{ A matrix containing the potential covariates to include in the model with one covariate per column and one individual per row}
  \item{get.OR.value}{ TRUE/FALSE: whether to return the OR values associated to the regression}
  \item{alpha}{ The alpha threshold to use for the OR confidence interval}
}
\details{
  This function will return results from the regression on each genomic region.

  If only two groups of individuals are present, a classical logistic regression is performed. 
  If more than two groups of individuals are present, a non-ordinal multinomial regression is performed, 
  comparing each group of individuals to the reference group indicated by the argument \code{ref.level}.

  The genetic score corresponding to CAST or WSS is directly computed in the function if \code{burden}="CAST" or "WSS". 
 
  If another genetic score is wanted, a matrix with one row per individual and one column per \code{genomic.region} 
  containing this score should be given to \code{burden}. In this situation, no bed matrix \code{x} is needed.
}
\value{
  A dataframe with one row per genomic region and at least two columns:
  \item{p.value}{ The p.value of the regression}
  \item{is.err}{ 0/1: whether there was a convergence problem with the regression}
  If \code{get.OR.value=TRUE}, additional columns are present:
  \item{OR}{ The OR value(s) associated to the regression. If there are more than two groups, there will be one OR value per group compared to the reference group}
  \item{l.lower}{ The lower bound of the confidence interval of each OR}
  \item{l.upper}{ The upper bound of the confidence interval of each OR}
}
\seealso{ \code{\link{power}}, \code{\link{CAST}}, \code{\link{WSS}} }
\examples{
#Data simulation: GRR values according to the SKAT formula
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov, 
                              select.gene="R1", n.case.groups = 2,
                              GRR.multiplicative.factor=2)

#Simulation of one group of 1,000 controls and two groups of 500 cases,
#50% of deleterious, 5 genomic regions are simulated.
x <- random.bed.matrix(genes.maf=Kryukov,
                       size = c(1000, 500, 500), baseline = c(0.001, 0.001),
                       GRR.matrix.del = GRR.del, prop.del = 0.5,
                       same.variant = FALSE, fixed.variant.prop = TRUE,
                       genetic.model = "multiplicative", replicates=5)

#Keep only rare variant
x <- filter.rare.variants(x, filter = "any", maf.threshold = 0.01, 
                          min.nb.snps = 5)
  							
#Simulation of covariates
X <- data.frame(sex = c(sample(0:1, 1000, TRUE, c(0.2,0.8)), 
                      sample(0:1, 1000, TRUE, c(0.8,0.2))), 
                u = runif(1000))
  				  
#Regression on 3 groups with the score WSS and OR values without covariates
burden.mlogit(x, burden = "WSS", maf.threshold=0.01, ref.level="0", 
              get.OR.value=TRUE, alpha = 0.05)
              
#Regression on 3 groups with the score WSS and OR values with the covariate "sex"
burden.mlogit(x, burden = "WSS", maf.threshold=0.01, ref.level="0", 
              get.OR.value=TRUE, alpha = 0.05, formula = ~ sex, data = X)
}   
