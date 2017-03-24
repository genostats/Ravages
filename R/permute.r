
permute <- function(x, centre, region, STAT, target = 10, B.max = 1e5, which.snps = rep(TRUE, ncol(x))) {
  Obs <- STAT(x, centre, region, which.snps)
  n.reg <- nlevels(region)
  A <- B <- integer(n.reg)
  b <- 0;
  while(b < B.max) {
    b <- b+1
    Perm <- STAT(x, sample(centre), region, which.snps);
    w.greater <- (Perm > Obs)
    A <- A + (!is.na(w.greater) & w.greater)
    B <- B + !is.na(w.greater)  
    # ne considérer que les snps qui sont dans un groupe où A < target
    keep <- (A < target)
    which.snps <- which.snps & keep[region]
  }
  data.frame(Obs = Obs, A = A, B = B, p = A/B)
}
