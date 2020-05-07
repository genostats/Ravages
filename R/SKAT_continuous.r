SKAT.continuous <- function(x, NullObject, genomic.region = x@snps$genomic.region, 
                     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                     estimation.pvalue = "kurtosis", debug = FALSE){
  x@snps$weights <- weights
  x <- select.snps(x, maf <= maf.threshold & maf > 0)
  x@snps$genomic.region <- as.factor(genomic.region)
  
  pheno <- matrix(NullObject$pheno, ncol = 1)
  
  #Matrix P1
  P1 <- NullObject$P1
  
  #(Y - pi_hat)
  ymp = NullObject$ymp
  
  #P-value for all regions
  res.allregions <- do.call(rbind, lapply(levels(x@snps$genomic.region), function(reg) get.parameters.pvalue.continuous(x, region = reg, P1 = P1, ymp = ymp, estimation.pvalue = estimation.pvalue)))
  res.final <- as.data.frame(res.allregions)
  colnames(res.final) <- colnames(res.allregions)
  rownames(res.final) <- levels(genomic.region)
  if(debug)
    res.final
  else
    res.final[, c("stat", "p.value")]
} 
  

#P-value by genomic region
get.parameters.pvalue.continuous <- function(x, region, P1, ymp, estimation.pvalue){
  x.genomic.region <- select.snps(x, genomic.region == region)
  G <- gaston::as.matrix(x.genomic.region) %*% diag(x.genomic.region@snps$weights)
  
  #Stat de test
  Q <- as.vector(ymp) %*% (G %*% t(G)) %*% as.vector(ymp)

  #Moments
  M <- .Call("moments", PACKAGE = "Ravages", G, P1)

  #P-valeur
  pval <- p.valeur.moments.liu(Q = Q, mu = M["mu"], sigma = M["sigma"], skewness = M["skewness"], kurtosis = M["kurtosis"], estimation.pvalue = estimation.pvalue)

  return(c(stat = Q, p.value = pval, mean = as.numeric(M["mu"]), M["sigma"], M["skewness"], M["kurtosis"]+3))
}


  
