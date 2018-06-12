# ----------- nvelle fonction
# pop.maf = un vecteur de mafs pour la population (m)
# OR = une matrice d'OR (n x m)  (on fait des efforts pour accepter vecteurs de longueur m (si n=1) ou n, ou matrices n x m')
# baseline = la prévalence de chacun des groupes (n)
# résultat = une matrice (n+1) x m (première ligne pour les témoins)
group.mafs <- function(pop.maf, OR, baseline) {
  if(!is.matrix(OR)) {
    if(length(OR) == length(pop.maf) & length(baseline) == 1)
      OR <- matrix(OR, nrow = 1)
    else if(length(OR) == length(baseline))
      OR <- matrix( rep_len(OR, length(pop.maf)*length(OR)), nrow = length(baseline))
    else
      stop("OR dimension mismatch")
  }

  if(nrow(OR) == length(baseline) & ncol(OR) < length(pop.maf)) 
    OR <- matrix( rep_len(OR, length(pop.maf)*nrow(OR)), nrow = nrow(OR))
  
  if(nrow(OR) != length(baseline) | ncol(OR) != length(pop.maf)) 
    stop("OR dimensions mismatch")

  if(length(baseline) != nrow(OR)) 
    baseline <- rep_len(baseline, nrow(OR))

  if(any(baseline < 0) | sum(baseline) > 1)
    stop("Unappropriate baseline values")

  p.t <- numeric(ncol(OR))
  for(i in 1:ncol(OR)) 
    p.t[i] <- p.tem(pop.maf[i], OR[,i], baseline)

  odds.t <- p.t/(1-p.t)
  F <- p.t
  for(i in 1:nrow(OR)) {
    odds.c <- OR[i,]*odds.t
    F <- rbind(F, odds.c/(1+odds.c))
  }

  if(nrow(F) == 2) 
    rownames(F) <- c("controls","cases")
  else 
    rownames(F) <- c("controls", sprintf("cases_%d", 1:(nrow(F)-1)))
 
  F
}


# pour calculer la maf témoins
# OR = un vecteur d'OR [un pour chaque groupe]
p.tem <- function(p, OR, baseline) {
  if(all(baseline == 0)) return(p); # si maladie très rare, maf = pop maf
  if(length(OR) == 1) { # eq du second degré
    if(OR == 1) return(p)
    alpha <- (OR-1)*(1-baseline)
    beta  <- 1+(OR-1)*(baseline - p)
    return((sqrt(beta**2 + 4*alpha*p) - beta)/(2*alpha))
  }
  B <- sum(baseline)
  f <- function(pt) (1-B + sum( baseline*OR/(1+(OR-1)*pt)))*pt - p
  return( uniroot(f, c(0,1), tol=1e-8)$root )
}
