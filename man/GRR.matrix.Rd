\name{GRR.matrix}
\alias{GRR.matrix}
\title{
  GRR matrix for genetic data simulation
}
\description{
  Computes a GRR matrix based on parameters for variants and genetic effects
}
\usage{
GRR.matrix(genes.maf = Kryukov, n.case.groups = 2, 
           GRR = c("SKAT", "constant", "variable"), 
           GRR.value, GRR.function, 
           GRR.multiplicative.factor = 2, select.gene)
}
\arguments{
  \item{genes.maf}{ A data frame containing the MAF in the general population (column \code{maf}) for variants grouped in genes (column \code{gene}). By default, \code{maf} from the file Kryukov are used}
  \item{n.case.groups}{ The number of cases groups (set as 2 by default), i.e. the number of groups where variants will have a GRR greater than 1}
  \item{GRR}{ How to calculate the GRR }
  \item{GRR.value}{ GRR value if \code{GRR="constant"}}
  \item{GRR.function}{ A function indicating how to calculate the GRR depending on MAF in the general population, only needed if \code{GRR="variable"}}
  \item{GRR.multiplicative.factor}{ A vector of size (\code{n.case.groups}-1) containing the multiplicative factor for the GRR for each group of cases compared to the first group of cases}
  \item{select.gene}{ The name of the gene(s) to select in the file \code{genes.maf} if multiple genes are present. If missing, the first level of \code{genes.maf$gene} is kept.}
}
\value{
  A matrix containing the GRR for each variant (in columns) in each group of cases (in lines)
}
\details{
  The GRR can be computed in three ways using the argument \code{GRR}.

  If \code{GRR="constant"}, the same GRR is given to all the variants with the value given to the argument \code{GRR.value}.
  If \code{GRR="SKAT"}, the GRR are calculating using the formula from the paper presenting the SKAT method and thus depend on MAF.
  If \code{GRR="variable"}, the GRR are calculating using a function given by the user to \code{GRR.function} depending only on the MAF in the general population.  

  The argument \code{multiplicative.factor} contain \code{n.case.groups}-1 values; if \code{multiplicative.factor=1}, GRR will be the same between the different groups of cases.  

  The data frame Kryukov (used by default) and GnomADgenes (containing MAF in the NFE population) can be used as \code{genes.maf}.

  \code{GRR.matrix} returns a matrix that can be used in other simulation functions such as \code{random.bed.matrix}.
}
\seealso{ \code{\link{random.bed.matrix}}, \code{\link{power}} }
\examples{
#GRR calculated on the MAF from the first unit of the file Kryukov
#using the formula from the SKAT paper, with the second group of cases 
#having GRR values twice as high as the first one
GRR.del <- GRR.matrix(GRR = "SKAT", genes.maf = Kryukov,
                      GRR.multiplicative.factor=2, 
                      select.gene = "R1")
}
