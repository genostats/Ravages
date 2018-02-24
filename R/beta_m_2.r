Beta.M <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
                   which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e6) {
  r <- .Call('oz_beta_m', PACKAGE = 'oz', x@bed, which.snps, as.factor(genomic.region), as.factor(group), target, B.max) 
  res <- data.frame( genomic.region = levels(genomic.region), r )
  return(res)
}

Beta.M.rect <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
                   which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e6) {
  r <- .Call('oz_beta_m', PACKAGE = 'oz', x@bed, which.snps, x@p, as.factor(genomic.region), as.factor(group), target, B.max) 
  res <- data.frame( genomic.region = levels(genomic.region), r )
  return(res)
}

Beta.M.exact <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
                   which.snps = rep(TRUE, ncol(x)), regions.to.test = levels(genomic.region)) {
  if(!is.factor(genomic.region))
    genomic.region <- as.factor(genomic.region)

  g <- match(regions.to.test, levels(genomic.region))
  if(any(is.na(g))) 
    stop("Unknown region")
  r <- .Call('oz_ex_beta_m', PACKAGE = 'oz', x@bed, which.snps, genomic.region, as.factor(group), g) 
  data.frame( genomic.region = regions.to.test, r )
}
