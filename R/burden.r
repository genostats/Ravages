burden <- function(x, NullObject, genomic.region = x@snps$genomic.region, burden, maf.threshold = 0.5, get.effect.size = FALSE, alpha = 0.05, cores = 10, verbose = TRUE){
  #Test if NullObject de bon type
  if("P1" %in% names(NullObject)) stop("'NullObject' has been generated with wrong 'RVAT' in 'NullObject.parameters()'") 
 
  if(missing(x)) x <- NULL
  if(NullObject$pheno.type == "categorial"){
    if(verbose) cat("Categorial phenotype \n")
    res <- burden.mlogit(x = x, NullObject = NullObject, genomic.region = genomic.region, burden = burden, maf.threshold = maf.threshold, get.effect.size = get.effect.size, alpha = alpha, cores = cores)
  }
  if(NullObject$pheno.type == "continuous"){
    if(verbose) cat("Continuous phenotype \n")
    res <- burden.continuous(x = x, NullObject = NullObject, genomic.region = genomic.region, burden = burden, maf.threshold = maf.threshold, get.effect.size = get.effect.size, alpha = alpha, cores = cores)
  }
  return(res)
}
