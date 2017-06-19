
permute <- function(x, centre, region, STAT, target = 10, B.max = 1e5, which.snps = rep(TRUE, ncol(x)), pearson=FALSE) {
  # mettre les niveaux du facteur 'region' dans l'ordre où ils apparaissent
  region1 <- as.character(region)
  region <- factor(region1,  unique(region1))
  # pas de NA dans which.snps, ça peut perturber STAT
  which.snps[ is.na(which.snps) ] <- FALSE
  # mettre à FALSE les snps dont le niveau est à NA
  which.snps[ is.na(region) ] <- FALSE

  Obs <- STAT(x, centre, region, which.snps)
  n.reg <- nlevels(region)
  A <- B <- S <- SC <- SC3 <- SC4 <- integer(n.reg)
  b <- 0;
  
  #Calcul de deux statistiques permutées si on veut réaliser une approximation par la loi de Pearson
  if(pearson==TRUE){
    a1 <- STAT(x, sample(centre), region, which.snps)
    a2 <- STAT(x, sample(centre), region, which.snps)
    K <- ( a1+a2 ) / 2
    L <- abs(a1-a2)
    L[which(L==0)] <- 1
    Obs <- (Obs - K) / L
  }
              
  while(b < B.max) {
    b <- b+1
    ifelse(pearson==TRUE, Perm <- ( (STAT(x, sample(centre), region, which.snps) - K ) / L) , Perm <- STAT(x, sample(centre), region, which.snps) );
    w.greater <- (Perm >= Obs)
    A <- A + (!is.na(w.greater) & w.greater)
    B <- B + !is.na(w.greater)  
    # ne considérer que les snps qui sont dans un groupe où A < target
    keep <- (A < target)
    which.snps <- which.snps & keep[region]
    Perm[ is.na(Perm) ] <- FALSE
    S <- S + Perm
    SC <- SC + Perm^2
    SC3 <- SC3 + Perm^3
    SC4 <- SC4 + Perm^4
  }
  #Calcul moyenne et variance
  M <- S/B
  V <- (1/ (B-1) ) * (SC - (S^2)/B)
  
  #Calcul des moments
  M2 <- (1/B) * (SC - (S^2)/B)
  M3 <- (1/B) * SC3 - 3 * M2 * M - M^3
  M4 <- (1/B) * SC4 - 4 * M * (1/B) * SC3 + 6 * (1/B) * SC * M^2 - 3 * M^4
  Sk <- M3 / (M2^(3/2))
  Ku <- M4 /(M2^2)
  
  p <- ifelse(B==B.max, pnorm( ((Obs - M) / sqrt(V)), lower.tail=FALSE), A/B)

  #Calcul de la p valeur avec approximation de Pearson
  if(pearson==TRUE){
    p <- ifelse(B==B.max, ppearson(Obs, pearsonFitM(M, V, Sk, Ku), lower.tail=FALSE), A/B)
  }
  
  pval <- data.frame(Obs = Obs, A = A, B = B, Moyenne = M, Variance = V, M2, M3, M4, Sk, Ku, p = p)
  
  #Calcul de la p valeur de pearson quand NA
  for (i in which(pval$p=="NaN")){
    pval[i,"p"] <- ppearson(pval$Obs[i], pearsonFitM(pval$Moyenne[i], pval$Variance[i], pval$Sk[i], pval$Ku[i]), lower.tail=FALSE) 
    }
        
  return(pval)
}
