SKAT <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, Pi, weights = x@p**(-24), maf.threshold = 0.01) {

  which.snps <- (x@p < maf.threshold) & (x@p > 0)
  genomic.region <- as.factor(genomic.region)
  group <- as.factor(group)

  # matrice des proba d'appartenir au group               
  if(missing(Pi)) {
    a <- table(group)/nrow(x)
    Pi <- matrix( a, ncol = nlevels(group), nrow = nrow(x), byrow = TRUE)
  }
  Pi <<- Pi
  B <- .Call('skat', PACKAGE = "Ravages", x@bed, which.snps, genomic.region, group, Pi, weights, 1, 0);
  return(B);
}

