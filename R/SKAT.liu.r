SKAT.liu <- function(x, NullObject, genomic.region = x@snps$genomic.region, 
                     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                     useskew = TRUE){
  x@snps$weights <- weights
  x <- select.snps(x, maf <= maf.threshold & maf > 0)
  genomic.region <- as.factor(genomic.region)
  
  group <- NullObject$group
  
  lev <- levels(group) 
  #Number of groups and individuals
  ngpe <- nlevels(group) ; nind <- length(group)
  #Number of individuals by group
  n <- as.vector(table(group))
  
  #Matrix pi
  Pi <- NullObject$Pi.data
  #Matrix X
  X <- NullObject$X
  
  #Matrix of Pi in all groups expect first one (by default)
  P0 <- P.mat( matrix(Pi[, -1], ncol=ngpe-1), X)
  #Matrix of Pi in all groups
  P1 <- P.mat2(Pi, X)

  #Matrix of (YY - Pi^) with YY indicatrice variable in each group 
  YY <- sapply(lev, function(l) as.numeric(group == l))
  ymp = YY -  Pi
  
  #P-value for all regions
  res.allregions <- t(sapply(levels(genomic.region), function(reg) get.pvalue.genomic.region(x, region = reg, P1 = P1, ymp = ymp, n = n, useskew = useskew)))
  return(res.allregions)
} 
  

#P-value by genomic region
get.pvalue.genomic.region <- function(x, region, P1, ymp, n, useskew = TRUE){
  ngpe <- length(n)
  x.genomic.region <- select.snps(x, genomic.region == region)
  #Matrix of weighted genotypes
  G <- gaston::as.matrix(x.genomic.region) %*% diag(x.genomic.region@snps$weights)
  #Kernel matrix
  K <- G %*% t(G)
  K.L <- lapply(1:ngpe, function(gpe) (1/n[gpe])*K)
  #Matrix by bloc
  K.bloc <- block.diag(K.L)
  
  #Stat de test
  Q <- as.vector(ymp) %*% K.bloc %*% as.vector(ymp)
  
  #Moments
  #M <- moments1(A = K.bloc, P = P1)
  M <- .Call("moments", PACKAGE = "Ravages", K.bloc, P1)
  
  #P-valeur
  pval <- p.valeur(Q = Q, moments = M, useskew = useskew)
  
  return(c(stat = Q, p.value = pval))
}
  
  
