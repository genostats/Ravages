SKAT <- function(x, NullObject, genomic.region = x@snps$genomic.region,
                 weights = (1 - x@snps$maf)**24, maf.threshold = 0.5, 
                 method = "size.based", params.bootstrap){
  if(method == "size.based"){
    if(nrow(x)<2000) method = "bootstrap"
    else  method = "liu"
  }else{
    if(!(method %in% c("bootstrap", "liu", "liu.kurtosis"))) stop("Wrong 'method' specified")
  }

  if(method == "bootstrap"){
    if(is.missing(params.bootstrap))
      params.bootstrap <- list(perm.target = 100, perm.max=1e4, debug = FALSE)
    res <- SKAT.bootstrap(x, NullObject, genomic.region = genomic.region,
                          weights = weights, maf.threshold = maf.threshold,
                          perm.target = params$perm.target, 
                          perm.max = params$perm.max, debug = params$debug)
  }

  if(method == "liu")
    res <- SKAT.Liu(x, NullObject, genomic.region = genomic.region, 
                    weights = weights, maf.threshold = maf.threshold,
                    useskew = TRUE)

  if(method == "liu.kurtosis")
    res <- SKAT.Liu(x, NullObject, genomic.region = genomic.region,
                    weights = weights, maf.threshold = maf.threshold,
                    useskew = FALSE)

  return(list(method = method, results = res))
}
