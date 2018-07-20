\name{power.burden}
\alias{power.burden}
\title{
  Simulation and power calculation
}
\description{
  Performs simulations and calculate the power of the logistic/multinomial regression on a genetic score as CAST or WSS
}
\usage{
  power.burden(alpha = 2.5e-06, filter = c("whole", "controls", "any"), min.nb.snps = NULL,
  			   maf.threshold = 0.01, CAST = TRUE, WSS = TRUE, other.score=NULL,
  			   pooled.analysis = TRUE, non.pooled.analysis = TRUE, analysis.by.group = TRUE, file.pop.maf = Kryukov,
  			   size = c(1000, 500, 500), baseline = c(0.001, 0.001), replicates = 1000, select.gene = NULL, 
  			   same.variant = FALSE, genetic.model = c("multiplicative", "general", "recessive", "dominant"),
  			   GRR.matrix, GRR.matrix.pro = NULL, prop.del = 0.5, prop.pro = 0, reflevel = "0")
}
\arguments{
  \item{alpha}{ The significance threshold to calculate the power}
  \item{filter}{ How to filter rare variants (see function \code{filter.rare.variants}}
  \item{maf.threshold}{ The MAF threshold used for the rare variant definition}
  \item{min.nb.snps}{ The minimum number of variants needed in a genomic region to keep it}
  \item{CAST}{ TRUE/FALSE: whether to calculate the power using the CAST score}
  \item{WSS}{ TRUE/FALSE: whether to calculate the power using the WSS score}
  \item{other.score}{ A function to calculate the genetic score, depending only on the bed matrix x simulated in \code{power.burden} if another score than CAST or WSS is wanted}
  \item{pooled.analysis}{ TRUE/FALSE: whether to calculate the power of the pooled analysis}
  \item{non.pooled.analysis}{ TRUE/FALSE: whether to calculate the power of the non-pooled analysis}
  \item{analysis.by.group}{ TRUE/FALSE: whether to calculate the power of the analysis by group}
  \item{file.pop.maf}{ A file containing the MAF in the general population (column maf) for variants within a gene (column gene), by default the file Kryukov from the package is used}
  \item{size}{ A vector giving the size of each group of individuals. By default: one group of 1,000 controls and two groups of 500 cases}
  \item{baseline}{ A vector giving the baseline of each group of cases, by default 0.001 for the two groups of cases}
  \item{replicates}{ The number of simulations to perform for power calculation, 1,000 by default}
  \item{select.gene}{ Which gene to choose in the column "gene" from \code{file.pop.maf} if multiple genes are present}
  \item{same.variant}{ TRUE/FALSE: whether the causal variants are the same in the different groups of cases}
  \item{genetic.model}{ The genetic model of the disease}
  \item{GRR.matrix}{ A list containing the GRR matrix associated to the heterozygous compared to the homozygous for the reference allele for deleterious variants. An additional GRR matrix associated to the homozygous for the alternate allele is needed if \code{genetic.model}="general"}
  \item{GRR.matrix.pro}{ The same argument as \code{GRR.matrix} but for protective variants}
  \item{prop.del}{ The proportion of deleterious variants in the cases, 0.5 by default}
  \item{prop.pro}{ The proportion of protective variants in the cases, 0 by default}
  \item{reflevel}{ The reference group in the regression, considered as the controls group by default}
}
\details{
  This function performs simulations before calculating the power associated to the regression models. It depends on multiple functions from the package.
  The regression model relies on a genetic score. The function can directly compute the CAST or WSS score if \code{CAST}=TRUE or \code{WSS}=TRUE. Another score can be specified in \code{other.score} which should be a function calculating the genetic score based solely on the bed matrix simulated in the programm.
  Three types of analysis are possible: (i) a pooled analysis where all the cases are compared as one group to the controls using a logistic regression ; (ii) a non pooled analysis where the different groups of cases are considered and compared to the grouo \code{reflevel} and a non-ordinal multinomial regression model is used ; (iii) an analysis by group where each group of cases is compared to the controls using a logistic regression.
  The MAF in the general population can be definied using the files Kryukov or GnomADgenes both available with the package.
  The genetic model of the disease needs to be specified in this function.
  If \code{model}="general", there is no link between the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele. Therefore, the user has to give two matri$
  If \code{model}="multiplicative", we assume that the the GRR associated to the homozygous for the alternate allele is the square of the GRR associated to the heterozygous: only one GRR matrix is needed.
  If \code{model}="dominant", we assume that the GRR associated to the heterozygous and the GRR associated to the homozygous for the alternate allele are equal: only one GRR matrix is needed.
  If \code{model}="recessive", we assume that the GRR associated to the heterozygous is equal to 1: the GRR given is the one associated to the homozygous for the alternate allele.
  The function compute.GRR.matrix from this package can help to get the GRR matrices which have to be passed as argument \code{GRR.matrix} and \code{GRR.matrix.pro}. These two matrices contain GRR value considering all the variants causal, these values are then attributed to the observed sampled causal variants in the simulations.
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
  GRR.del <- compute.GRR.matrix(GRR = "SKAT", file.pop.maf = Kryukov, n.case.groups = 2, GRR.multiplicative.factor=2, select.gene = "R1")
                                  
  #Simulation of one group of 1,000 controls and two groups of 500 cases, each one with a prevalence of 0.001.
  #50% of deleterious variants having a GRR calculated with the SKAT formula based on MAF in the file Kryukov are present and they are different between the two groups of cases.
  #100 genomic regions are simulated to calculate the power of CAST and WSS at a significance threshold of 2.5e-6.
  #Only non-monomorphic variants with a MAF lower than 1% in any of the groups are kept. Only regions with at least 5 variants are kept.
  power.burden(alpha = 0.001, WSS = TRUE, CAST = TRUE, maf.threshold = 0.01, filter = "any", min.nb.snps = 5, pooled.analysis = TRUE, non.pooled.analysis = TRUE, analysis.by.group = TRUE,
               file.pop.maf = Kryukov, select.gene = "R1", size = c(1000,500,500), baseline = c(0.001,0.001), same.variant = FALSE, GRR.matrix = GRR.del,
               genetic.model = "multiplicative", prop.del = 0.5, prop.pro = 0, reflevel = "0", replicates = 100)  
  }
