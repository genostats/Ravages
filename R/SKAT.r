SKAT <- function(x, NullObject, genomic.region = x@snps$genomic.region,
                 weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                 get.moments = "size.based", estimation.pvalue = "kurtosis", 
                 params.sampling){
  if(get.moments == "size.based"){
    if(nrow(x)<2000) get.moments = "bootstrap"
    else  get.moments = "theoretical"
  }else{
    if(!(method %in% c("permutations", "bootstrap", "theoretical"))) stop("Wrong 'method' specified")
  }

  if(method == "bootstrap"){
    if(missing(params.sampling))
      params.sampling <- list(perm.target = 100, perm.max=1e4, debug = FALSE)
    res <- SKAT.bootstrap(x, NullObject, genomic.region = genomic.region,
                          weights = weights, maf.threshold = maf.threshold,
                          perm.target = params.bootstrap$perm.target, 
                          perm.max = params.bootstrap$perm.max, 
                          debug = params.bootstrap$debug, estimation.pvalue = estimation.pvalue)
  }
  
  if(method == "permutations"){
    if(missing(params.sampling))
      params.sampling <- list(perm.target = 100, perm.max=1e4, debug = FALSE)
    res <- SKAT.permutations(x, NullObject, genomic.region = genomic.region,
                             weights = weights, maf.threshold = maf.threshold,
                             perm.target = params.bootstrap$perm.target, 
                             perm.max = params.bootstrap$perm.max, 
                             debug = params.bootstrap$debug, estimation.pvalue = estimation.pvalue)
  }

  if(method == "theoretical")
    res <- SKAT.Liu(x, NullObject, genomic.region = genomic.region, 
                    weights = weights, maf.threshold = maf.threshold,
                    useskew = TRUE)

  return(list(method = method, results = res))
}
