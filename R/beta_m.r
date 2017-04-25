Beta.M <- function(x, centre, region, which.snps = rep(TRUE, ncol(x))) {
  L <- alleles.by.group(x, centre, which.snps)
  Betam <- colSumsSq(L$minor)/colSums(L$minor) + colSumsSq(L$major)/colSums(L$major)
  R <- sum_by_index(Betam, region[which.snps])
  R[ R == 0 ] <- NA # les groupes où which.snps = FALSE pour tous les SNPs !
  R
}
