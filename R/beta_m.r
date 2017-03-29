Beta.M <- function(x, centre, region, which.snps = rep(TRUE, ncol(x))) {
  L <- alleles.by.group(x, centre, which.snps)
  Betam <- colSumsSq(L$minor)/colSums(L$minor) + colSumsSq(L$major)/colSums(L$major)
  .Call("oz_sum_by_group", PACKAGE = "oz", Betam, region[which.snps])
}
