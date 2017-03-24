
minor.alleles.by.group <- function(x, group, which.snps = rep(TRUE, ncol(x))) {
  if(length(which.snps) != ncol(x) | typeof(which.snps) != "logical")
    stop("which.snps must be a logical vector of length = ncol(x)")

  if(length(group) != nrow(x) | !is.factor(group))
    stop("group must be a factor of length = nrow(x)")

  inverse <- (x@p > 0.5)
 
  C <- .Call('oz_alt_alleles_by_factor', PACKAGE = 'oz', x@bed, which.snps, group, inverse)
  rownames(C) <- levels(group)
  id <- x@snps$id[which.snps]
  if(!anyDuplicated(id)) colnames(C) <- id
  C
}

major.alleles.by.group <- function(x, group, which.snps = rep(TRUE, ncol(x))) {
  if(length(which.snps) != ncol(x) | typeof(which.snps) != "logical")
    stop("which.snps must be a logical vector of length = ncol(x)")

  if(length(group) != nrow(x) | !is.factor(group))
    stop("group must be a factor of length = nrow(x)")

  inverse <- (x@p <= 0.5)
 
  C <- .Call('oz_alt_alleles_by_factor', PACKAGE = 'oz', x@bed, which.snps, group, inverse)
  rownames(C) <- levels(group)
  id <- x@snps$id[which.snps]
  if(!anyDuplicated(id)) colnames(C) <- id
  C
}
