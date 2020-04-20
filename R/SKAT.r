SKAT <- function(x, NullObject, genomic.region = x@snps$genomic.region,
                 weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                 get.moments = "size.based", estimation.pvalue = "kurtosis", 
                 params.sampling, cores, debug = FALSE){
  if(get.moments == "size.based"){
    if(nrow(x)<2000) get.moments = "bootstrap"
    else  get.moments = "theoretical"
  }else{
    if(!(get.moments %in% c("permutations", "bootstrap", "theoretical"))) stop("Wrong 'get.moments' specified")
  }

  if(get.moments == "bootstrap"){
    if(missing(params.sampling))
      params.sampling <- list(perm.target = 100, perm.max=5e4)
    res <- SKAT.bootstrap(x, NullObject, genomic.region = genomic.region,
                          weights = weights, maf.threshold = maf.threshold,
                          perm.target = params.sampling$perm.target, 
                          perm.max = params.sampling$perm.max, 
                          debug = debug, estimation.pvalue = estimation.pvalue)
  }
  
  if(get.moments == "permutations"){
    if(missing(params.sampling))
      params.sampling <- list(perm.target = 100, perm.max=5e4)
    res <- SKAT.permutations(x, NullObject, genomic.region = genomic.region,
                             weights = weights, maf.threshold = maf.threshold,
                             perm.target = params.sampling$perm.target, 
                             perm.max = params.sampling$perm.max, 
                             debug = debug, estimation.pvalue = estimation.pvalue)
  }

  if(get.moments == "theoretical"){
    if(missing(cores)) cores <- 10
    res <- SKAT.theoretical(x, NullObject, genomic.region = genomic.region, 
                            weights = weights, maf.threshold = maf.threshold,
                            estimation.pvalue = estimation.pvalue, cores = cores, debug = debug)
  }

  return(list(get.moments = get.moments, results = res))
}
