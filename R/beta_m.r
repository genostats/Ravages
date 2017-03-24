Beta.M <- function(x, centre, region, which.snps = rep(TRUE, ncol(x))) {
  A <- minor.alleles.by.group(x, centre, which.snps)
  B <- major.alleles.by.group(x, centre, which.snps)
  Betam <- colSums(A**2)/colSums(A) + colSums(B**2)/colSums(B)
  .Call("oz_sum_by_group", PACKAGE = "oz", Betam, region[which.snps])
}
