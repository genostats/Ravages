C.ALPHA <- function(x, centre, region, which.snps = rep(TRUE, ncol(x))) {
  Nc <- minor.alleles.by.group(x, centre, which.snps)
  N <- colSums(Nc)
  Mc <- major.alleles.by.group(x, centre, which.snps)
  M <- colSums(Mc)
  p <- sweep( Nc + Mc, 2, (N + M ), "/" )
  Sc <- (sweep(p, 2, N, "*") - Nc)**2 - sweep(p*(1-p), 2, N, "*")
  Ca <- colSums(Sc)
  sp2 <- colSums(p^2)
  Vc <- 2 * N * ( (N-3) * sp2^2 + (N-1) * sp2 - 2 * (N-2) * colSums(p^3) )
  Ca_Sum <- .Call("oz_sum_by_group", PACKAGE = "oz", Ca, region[which.snps])
  Vc_Sum <- .Call("oz_sum_by_group", PACKAGE = "oz", Vc, region[which.snps])
  Ca_Sum/sqrt(Vc_Sum)
}