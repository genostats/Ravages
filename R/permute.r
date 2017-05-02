
permute <- function(x, centre, region, STAT, target = 10, B.max = 1e5, which.snps = rep(TRUE, ncol(x))) {
  # mettre les niveaux du facteur 'region' dans l'ordre où ils apparaissent
  region1 <- as.character(region)
  region <- factor(region1,  unique(region1))
  # pas de NA dans which.snps, ça peut perturber STAT
  which.snps[ is.na(which.snps) ] <- FALSE
  # mettre à FALSE les snps dont le niveau est à NA
  which.snps[ is.na(region) ] <- FALSE

  Obs <- STAT(x, centre, region, which.snps)
  n.reg <- nlevels(region)
  A <- B <- S <- SC <- integer(n.reg)
  b <- 0;
  while(b < B.max) {
    b <- b+1
    Perm <- STAT(x, sample(centre), region, which.snps);
    w.greater <- (Perm >= Obs)
    A <- A + (!is.na(w.greater) & w.greater)
    B <- B + !is.na(w.greater)  
    # ne considérer que les snps qui sont dans un groupe où A < target
    keep <- (A < target)
    which.snps <- which.snps & keep[region]
    Perm[ is.na(Perm) ] <- FALSE
    S <- S + Perm
    SC <- SC + Perm^2
  }
  #Calcul moyenne et variance
  M <- S/B
  V <- (1/ (B-1) ) * (SC - (S^2)/B)
  
  
  p <- ifelse(B==B.max, pnorm( ((Obs - M) / sqrt(V)), lower.tail=FALSE), A/B)
  
  data.frame(Obs = Obs, A = A, B = B, Moyenne = M, S=S, Sc=SC, Variance = V, p = p)
}
