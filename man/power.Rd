\name{power}
\alias{power}
\title{
  Simulation and power calculation
}
\description{
  Performs simulations and calculates the power of the logistic/multinomial regression on a genetic score as CAST or WSS
}
\usage{
power(alpha = 2.5e-06, filter = c("whole", "controls", "any"), 
      min.nb.snps, maf.threshold = 0.01, 
      CAST = c(TRUE, FALSE), WSS = c(TRUE, FALSE), other.score,
      pooled.analysis = c(TRUE, FALSE), 
      non.pooled.analysis = c(TRUE, FALSE), 
      analysis.by.group = c(TRUE, FALSE), 
      model.pars=list(), 
      formula=NULL, data=NULL, ref.level = "0",
      model.pars=list(genes.maf=Kryukov, select.gene="R1", 
                      size=c(1000, 500, 500), baseline=c(0.001, 0.001), 
                      n.case.groups=2, replicates=1000,
                      GRR="SKAT", GRR.multiplicative.factor=2,
                      same.variant=FALSE, fixed.variant.prop = TRUE,
                      genetic.model="multiplicative",
                      prop.del=0.5, prop.pro=0))
}
\arguments{
  \item{alpha}{ The significance threshold, 2.5e-06 by default}
  \item{filter}{ On which group apply the filter rare variants (see function \code{filter.rare.variants})}
  \item{maf.threshold}{ The MAF threshold used for the rare variant definition}
  \item{min.nb.snps}{ The minimum number of variants needed to keep a genomic region}
  \item{CAST}{ TRUE/FALSE: whether to calculate the power using the CAST score}
  \item{WSS}{ TRUE/FALSE: whether to calculate the power using the WSS score}
  \item{other.score}{ A function to calculate the genetic score, depending only on the bed matrix simulated in \code{power} if another score than CAST or WSS is wanted}
  \item{pooled.analysis}{ TRUE/FALSE: whether to calculate the power of the pooled analysis, i.e. comparing all the cases to the controls}
  \item{non.pooled.analysis}{ TRUE/FALSE: whether to calculate the power of the non-pooled analysis, i.e. considering all the groups of individuals}
  \item{analysis.by.group}{ TRUE/FALSE: whether to calculate the power of the analysis by group, i.e. comparing each group of cases to the controls}
  \item{formula}{ An R formula corresponding to the regression model, if missing, only the genetic score will be included in the regression}
  \item{data}{ A matrix containing the potential covariates to include in the model with one covariate per column and one individual per row}
  \item{ref.level}{ The reference group in the regression, the controls group by default}
  \item{model.pars}{ A list containing the arugments from \code{random.bed.matrix}. The arguments from \code{GRR.matrix} or directly a matrix of GRR is also needed.}
}
\details{
  This function relies on multiple functions from \code{Ravages} and performs simulations before calculating the power associated to the regression models on a genetic score. 

  \code{power} can directly compute the CAST or WSS score if \code{CAST=TRUE} or \code{WSS=TRUE}. 
  Another score can be specified in \code{other.score} as a function calculating the genetic score based solely on the bed matrix simulated in the programm.  
  
  Three types of analysis are possible: (i) a pooled analysis where all the cases are compared as one group to the controls using a logistic regression ; 
  (ii) a non pooled analysis where the different groups of cases are considered using a non-ordinal multinomial regression model ; 
  (iii) an analysis by group where each group of cases is compared to the controls using a logistic regression.  
  
  To run the simulations, \code{power} needs to use the same arguments as \code{GRR.matrix} and \code{random.bed.matrix}:
  the arguments \code{genes.maf}, \code{select.gene}, \code{size}, \code{baseline}, \code{replicates}, \code{same.variant}, \code{fixed.variant.prop}, 
  \code{genetic.model}, \code{GRR.matrix.del}, \code{GRR.matrix.pro}, \code{prop.del}, \code{prop.pro} will be used by
  \code{power} and should therefore be given, otherwise the default values from \code{random.bed.matrix} will be used.
  The matrix of GRR values can be directly given to \code{GRR.matrix.del} (and \code{GRR.matrix.pro}) in the same format as in \code{random.bed.matrix},
  or this matrix can be computed directly in \code{power} by a call to \code{GRR.matrix}: in this situation, the arguments \code{n.case.groups}, \code{baseline},
  \code{GGR}, \code{GRR.value}/\code{GRR.function}, \code{GRR.mutiplicative.factor}, \code{select.gene} should be given, or the default values will be used.
  
  The MAF in the general population for the argument \code{model.pars$genes.maf} can be definied using the files \code{Kryukov} or \code{GnomADgenes} both available with the package.  
  
  The genetic model of the disease should be given as one parameter of \code{model.pars} to simulate data:
  If \code{genetic.model="general"}, there is no link between the GRR associated to the heterozygous genotype and the GRR associated to the homozygous alternate genotype. Therefore, the user has to give two matrices of GRR, one for each of these two genotypes.
  If \code{genetic.model="multiplicative"}, we assume that the the GRR associated to the homozygous alternated genotype is the square of the GRR associated to the heterozygous genotype: only one GRR matrix is needed.
  If \code{genetic.model="dominant"}, we assume that the GRR associated to the heterozygous genotype and the GRR associated to the homozygous alternate genotype are equal: only one GRR matrix is needed.
  If \code{genetic.model="recessive"}, we assume that the GRR associated to the heterozygous genotype is equal to 1: the GRR given is the one associated to the homozygous alternate genotype.  
}
\value{
  A dataframe with one row per analysis for each test (CAST, WSS or another score): pooled.analysis, each.cases.vs.tem, non.pooled.analysis; and three columns:
  \item{power}{ The power of the analysis}
  \item{se}{ The standard error of the power}
  \item{nb.replicates}{ The number of simulations without error in the regression model for the power estimation}
}
\seealso{ \code{\link{GRR.matrix}}, \code{\link{random.bed.matrix}}, \code{\link{CAST}}, \code{\link{WSS}} }
\examples{
#Simulation of one group of 1,000 controls and two groups of 500 cases, 
#50% of the variants are deleterious, 100 genomic regions are simulated
#GRR calculated using the SKAT formula
#Power is calculated at a significance threshold of 2.5e-6
power(alpha = 0.001, WSS = TRUE, CAST = TRUE, maf.threshold=0.01,
      filter = "any", min.nb.snps = 5,
      pooled.analysis = TRUE, non.pooled.analysis = TRUE, 
      analysis.by.group = TRUE, ref.level="0",
      formula=NULL, data=NULL,
      model.pars=list(genes.maf=Kryukov, select.gene="R1", size=c(1000, 500, 500), 
                      baseline=c(0.001, 0.001), n.case.groups=2, replicates=100, 
                      GRR="SKAT", GRR.multiplicative.factor=2,
                      same.variant=FALSE, fixed.variant.prop = TRUE, 
                      genetic.model="multiplicative",
                      prop.del=0.5, prop.pro=0) )
}
