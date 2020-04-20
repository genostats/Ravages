#P-value by genomic region
get.parameters.pvalue.theoretical <- function(x, region, P1, ymp, n, estimation.pvalue){
  ngpe <- length(n)

  x.genomic.region <- select.snps(x, genomic.region == region)
  G <- gaston::as.matrix(x.genomic.region) %*% diag(x.genomic.region@snps$weights)
  G.L <- lapply(1:ngpe, function(gpe) sqrt(1/n[gpe]) * G)
  G.bloc <- .Call("block_diag", PACKAGE = "Ravages", G.L)

  #Stat de test
  Q <- as.vector(ymp) %*% (G.bloc %*% t(G.bloc)) %*% as.vector(ymp)

  #Moments
  M <- .Call("moments", PACKAGE = "Ravages", G.bloc, P1)

  #P-valeur
  pval <- p.valeur.moments.liu(Q = Q, mu = M["mu"], sigma = M["sigma"], skewness = M["skewness"], kurtosis = M["kurtosis"], estimation.pvalue = estimation.pvalue)

  return(c(stat = Q, p.value = pval, mean = as.numeric(M["mu"]), M["sigma"], M["skewness"], M["kurtosis"]+3))
}

