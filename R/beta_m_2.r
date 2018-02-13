Beta.M <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
                   which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e6) {
  r <- .Call('oz_beta_m', PACKAGE = 'oz', x@bed, which.snps, as.factor(genomic.region), as.factor(group), target, B.max) 
  res <- data.frame( genomic.region = levels(genomic.region), r )
  return(res)
}

ex_Beta.M <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
                   which.snps = rep(TRUE, ncol(x)), groups) {
  r <- .Call('oz_ex_beta_m', PACKAGE = 'oz', x@bed, which.snps, as.factor(genomic.region), as.factor(group), groups) 
  r
}
