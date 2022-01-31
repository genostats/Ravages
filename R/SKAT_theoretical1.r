SKAT.theoretical1 <- function(x, NullObject, genomic.region = x@snps$genomic.region, 
                     weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                     estimation.pvalue = "kurtosis", cores = 10, debug = FALSE){


  #Check between number of individuals
  if(nrow(x) != length(NullObject$group)) stop("Different number of individuals in 'x' and 'NullObject'")

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)
  
  if(any(table(genomic.region)==1)) 
    stop("All 'genomic.region' sould contain at least 2 variants, please use 'filter.rare.variants()' to filter the bed matrix")

  which.snps <- (x@snps$maf <= maf.threshold) & (x@snps$maf > 0)
  
  group <- NullObject$group
  Pi <- NullObject$Pi.data
  P <- NullObject$P1

  stat <- .Call('rvg_skatStats', PACKAGE = "Ravages", x@bed, which.snps, genomic.region, group, x@p, Pi, P, weights);

  nb.ind.in.groups <- as.vector(table(group))

  M <- mclapply(1:nlevels(genomic.region) - 1L, 
                function(g) .Call('rvg_skatMoments', PACKAGE = "Ravages", 
                      x@bed, which.snps, genomic.region, group, x@p, Pi, P, weights, nb.ind.in.groups, g), mc.cores = cores);

  M <- as.data.frame(cbind(stat, t(simplify2array(M))))

  M$p.value <- mapply(p.valeur.moments.liu, Q = M$stat, mu = M$mean, sigma = M$sigma, skewness = M$skewness, kurtosis = M$kurtosis, 
                                    estimation.pvalue = estimation.pvalue)

  M <- M[,c("stat", "p.value", "mean", "sigma", "skewness", "kurtosis")]
  M$kurtosis <- 3 + M$kurtosis
  rownames(M) <- levels(genomic.region)
  if(debug)
    M
  else
    M[,1:2]
}
  

  
