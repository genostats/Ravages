SKAT.liu <- function(x, group=x@ped$pheno, genomic.region = x@snps$genomic.region, weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, ref.level, data = NULL, formula = NULL){
  x@snps$weights <- weights
  x <- select.snps(x, maf <= maf.threshold & maf > 0)
  genomic.region <- as.factor(genomic.region)
  if(!is.factor(group)) stop("group should be a factor")
  if(ref.level != levels(group)[1]) stop("ref.level should be the first level of group")
  lev <- levels(group) 
  #Number of groups and individuals
  ngpe <- nlevels(group) ; nind <- length(group)
  #Number of individuals by group
  n <- as.vector(table(group))
  
  #Matrix pi
  if(!is.null(data)){
    pi <- Ravages:::Pi.matrix(group, data=data, formula=formula, ref.level=ref.level)
    if(!is.null(formula)) X <- cbind(1, data[, all.vars(formula)])
    else X <- cbind(1, data)
  }
  if(is.null(data)){
    pi <- matrix(rep(n/nind, each=nind), ncol=ngpe)
    X <- matrix(1, ncol=1, nrow=nind)
  }
  
  #Matrix of pi in all groups expect ref.level
  P0 <- P.mat( matrix(pi[, -1], ncol=ngpe-1), X)
  #Matrix of pi in all groups
  P1 <- P.mat2(pi, X)

  #Matrix of (YY - pi^) with YY indicatrice variable in each group 
  YY <- sapply(lev, function(l) as.numeric(group == l))
  ymp = YY -  pi
  
  #P-value for all regions
  res.allregions <- sapply(levels(genomic.region), function(reg) get.pvalue.genomic.region(x, region = reg, P1 = P1, ymp = ymp, n = n))
  return(res.allregions)
} 
  

#P-value by genomic region
get.pvalue.genomic.region <- function(x, region, P1, ymp, n){
  ngpe <- length(n)
  x.genomic.region <- select.snps(x, genomic.region == region)
  #Matrix of weighted genotypes
  G <- as.matrix(x.genomic.region) %*% diag(x.genomic.region@snps$weights)
  #Kernel matrix
  K <- G %*% t(G)
  K.L <- lapply(1:ngpe, function(gpe) (1/n[gpe])*K)
  #Matrix by bloc
  K.bloc <- block.diag(K.L)
  
  #Stat de test
  Q <- as.vector(ymp) %*% K.bloc %*% as.vector(ymp)
  
  #P-valeur
  pval <- p.valeur(A = K.bloc, P = P1, Q = Q)
  
  return(c(stat = Q, p.value = pval))
}
  
  