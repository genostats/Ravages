\name{burden}
\alias{burden}

\title{ Burden }

\description{ Compute burdens .... }
\usage{ burden(x, regions, weights.1, weights.2 = weights.1) }

\arguments{
  \item{x}{ A \code{\link{bed.matrix}} }
}

\details{ bla bla
}

\value{
blu blu
}

% \seealso{  }

\examples{
# Load data
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)  # 361 SNPs

# define three regions / weight everything to 1
reg <- factor( c(rep("a", 100), rep("b",100), rep("c",161) ) )
wei.0 <- rep(0, 361)
wei.1 <- rep(1, 361)
wei.2 <- rep(4, 361)
B <- burden(x, reg, wei.0, wei.1, wei.2)
dim(B)
}

% \keyword{ Burden }
