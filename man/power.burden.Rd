\name{power.burden}
\alias{power.burden}
\title{
  Simulation and power calculation
}
\description{
  Performs simulations and calculate the power of the logistic/multinomial regression on a genetic score as CAST or WSS
}
\usage{
  power.burden(alpha = 2.5e-06, filter = c("whole", "controls", "any"), 
  			   min.nb.snps = NULL, maf.threshold = 0.01, 
  			   CAST = c(TRUE, FALSE), WSS = c(TRUE, FALSE), other.score=NULL,
  			   pooled.analysis = c(TRUE, FALSE), non.pooled.analysis = c(TRUE, FALSE), 
  			   analysis.by.group = c(TRUE, FALSE), 
  			   file.pop.maf = Kryukov, select.gene=NULL,
  			   size = c(1000, 500, 500), baseline = c(0.001, 0.001), replicates = 1000, 
  			   same.variant = c(FALSE, TRUE), fixed.variant.prop = c(TRUE, FALSE),
  			   model = c("multiplicative", "general", "recessive", "dominant"),
  			   GRR.matrix, GRR.matrix.pro = NULL, 
  			   prop.del = 0.5, prop.pro = 0, covariates=NULL, reflevel = "0")
}
\arguments{
  \item{alpha}{ The significance threshold, 2.5e-06 bu default}
  \item{filter}{ How to filter rare variants (see function \code{filter.rare.variants}}
  \item{maf.threshold}{ The MAF threshold used for the rare variant definition}
  \item{min.nb.snps}{ The minimum number of variants needed to keep a genomic region}
  \item{CAST}{ TRUE/FALSE: whether to calculate the power using the CAST score}
  \item{WSS}{ TRUE/FALSE: whether to calculate the power using the WSS score}
  \item{other.score}{ A function to calculate the genetic score, depending only on the bed matrix simulated in \code{power.burden} if another score than CAST or WSS is wanted}
  \item{pooled.analysis}{ TRUE/FALSE: whether to calculate the power of the pooled analysis, i.e. comparing all the cases to the controls}
  \item{non.pooled.analysis}{ TRUE/FALSE: whether to calculate the power of the non-pooled analysis, i.e. considering all the groups of individuals}
  \item{analysis.by.group}{ TRUE/FALSE: whether to calculate the power of the analysis by group, i.e. comparing each group of cases to the controls}
  \item{file.pop.maf}{ A file containing the MAF in the general population (column maf) for variants within genes (column gene), by default the file Kryukov is used}
  \item{select.gene}{ Which gene to choose in the column "gene" from \code{file.pop.maf} if multiple genes are present}
  \item{size}{ A vector giving the size of each group of individuals. By default: one group of 1,000 controls and two groups of 500 cases}
  \item{baseline}{ A vector giving the baseline of each group of cases, by default 0.001 for the two groups of cases}
  \item{replicates}{ The number of simulations to perform for power calculation. 1,000 by default}
  \item{same.variant}{ TRUE/FALSE: whether the causal variants are the same in the different groups of cases}
  \item{fixed.variant.prop}{ TRUE/FALSE: whether the arguments \code{prop.del} and \code{prop.pro} correspond to the final proportion of deleterious and protective variants (TRUE) or to the probability for each variant of being deleterious or protective (FALSE)}
  \item{model}{ The genetic model of the disease}
  \item{GRR.matrix}{ A list containing the deleterious GRR matrix associated to the heterozygous compared to the homozygous for the reference allele. An additional GRR matrix associated to the homozygous for the alternate allele is needed if \code{model}="general"}
  \item{GRR.matrix.pro}{ Needed if \code{prob.pro>0}: the same argument as \code{GRR.matrix} but with protective GRR values}
  \item{prop.del}{ The proportion of deleterious variants in cases, 0.5 by default}
  \item{prop.pro}{ The proportion of protective variants in cases, 0 by default}
  \item{reflevel}{ The reference group in the regression, considered as the controls group by default}
  \item{covariates}{ A matrix containing the potential covariates to include in the model with one covariate per column and one individual per row}
}
\details{
  This function relies on multiple functions from Ravages and performs simulations before calculating the power associated to the regression models. 
  The regression model relies on a genetic score. \code{power.burden} can directly compute the CAST or WSS score if \code{CAST}=TRUE or \code{WSS}=TRUE. Another score can be specified in \code{other.score} as a function calculating the genetic score based solely on the bed matrix simulated in the programm.
  Three types of analysis are possible: (i) a pooled analysis where all the cases are compared as one group to the controls using a logistic regression ; (ii) a non pooled analysis where the different groups of cases are considered and compared to the group \code{reflevel} using a non-ordinal multinomial regression model is used ; (iii) an analysis by group where each group of cases is compared to the controls using a logistic regression.
  The MAF in the general population can be definied using the files Kryukov or GnomADgenes both available with the package.
  The genetic model of the disease needs to be specified in this function.
  If \code{model}="general", there is no link between the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele. Therefore, the user has to give two matrices of GRR, one for the heterozygous, the other for the homozygous for the alternate allele.
  If \code{model}="multiplicative", we assume that the the GRR associated to the homozygous for the alternate allele is the square of the GRR associated to the heterozygous: only one GRR matrix is needed.
  If \code{model}="dominant", we assume that the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele are equal: only one GRR matrix is needed.
  If \code{model}="recessive", we assume that the GRR associated to the heterozygous is equal to 1: the GRR given is the one associated to the homozygous for the alternate allele.
  The function compute.GRR.matrix from this package can help to get the GRR matrices which have to be passed as argument \code{GRR.matrix} (and \code{GRR.matrix.pro}). These two matrices contain GRR value considering all the variants deleterious (and protective), these values are then attributed to the observed sampled deleterious variants (and protective) in the simulations.
}
\value{
  A dataframe with one row per analysis for each test (CAST, WSS or another score): pooled.analysis, each.cases.vs.tem, non.pooled.analysis; and three columns:
  \item{power}{ The power of the analysis}
  \item{se}{ The standard error of the power}
  \item{nb.replicates}{ The number of simulations without error in the regression model for the power estimation}
}
\author{
  Ozvan Bocher and Herve Perdy
}
\examples{
  #GRR matrix
  GRR.del <- compute.GRR.matrix(GRR = "SKAT", file.pop.maf = Kryukov, n.case.groups = 2, 
  								GRR.multiplicative.factor=2, select.gene = "R1")
                                 
  #Simulation of one group of 1,000 controls and two groups of 500 cases, 
  #each one with a prevalence of 0.001.
  #50% of the variants are deleterious with GRR calculated with the SKAT formula based on MAF
  #in the file Kryukov, and they are different between the two groups of cases.
  #100 genomic regions are simulated to calculate the power of CAST and WSS 
  #with a significance threshold of 2.5e-6.
  #Only non-monomorphic variants with a MAF lower than 1% in any of the groups 
  #and regions with at least 5 variants are kept.
  power.burden(alpha = 0.001, WSS = TRUE, CAST = TRUE, maf.threshold = 0.01, filter = "any", 
  			   min.nb.snps = 5, pooled.analysis = TRUE, non.pooled.analysis = TRUE, 
			   analysis.by.group = TRUE, file.pop.maf = Kryukov, select.gene = "R1", 
			   size = c(1000,500,500), baseline = c(0.001,0.001), 
			   same.variant = FALSE, fixed.variant.prop = TRUE, GRR.matrix = GRR.del,
               model = "multiplicative", prop.del = 0.5, prop.pro = 0, 
               reflevel = "0", replicates = 100)  
  }
