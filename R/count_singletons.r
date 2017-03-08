
count.singletons <- function(x) {
  w1 <- (x@snps$N1 == 1) & (x@snps$N2 == 0) 
  w2 <- (x@snps$N1 == 1) & (x@snps$N0 == 0) 
 
  C <- .Call('oz_count_alternative_alleles', PACKAGE = 'oz', x@bed, w1 | w2, w2 )
  names(C) <- x@ped$id
  C
}
