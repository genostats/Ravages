RAVA.FIRST <- function(x, variant.scores = NULL, ref.level, filter = c("whole", "controls", "any"), maf.threshold = 0.01, min.nb.snps = 2, min.cumulative.maf, group, cores = 10, burden = TRUE, H0.burden, burden.parameters, SKAT = TRUE, H0.SKAT, SKAT.parameters, verbose = TRUE){
  ###Attribute CADD regions to variants
  if(verbose) cat("Attribution of CADD regions\n")
  x <- set.CADDregions(x)
  ###Filtering based on CADD scores
  x.filter <- filter.adjustedCADD(x, variant.scores = variant.scores, ref.level = ref.level, filter = filter, maf.threshold = maf.threshold, min.nb.snps = min.nb.snps, min.cumulative.maf = min.cumulative.maf, group = group, verbose = verbose)
  if(verbose) cat("Filtering of rare variants within CADD regions\n")
  ###Association
  if(burden){
    if(missing(burden.parameters)){ burden.parameters <- list(burden.function = WSS, get.effect.size = F) }
    if(missing(H0.burden)) stop("'H0.burden should be given and can be obtained using 'NullObject.parameters()'")
    if(verbose) cat("Burden test\n")
    x.burden <- burden.subscores(x.filter, H0.burden, burden.function = burden.parameters$burden.function, maf.threshold = maf.threshold, get.effect.size = burden.parameters$get.effect.size, alpha = 0.05, cores = cores, verbose = verbose)
  }else{
    x.burden <- NULL
  }
  if(SKAT){
    if(missing(SKAT.parameters)){ SKAT.parameters <- list(weights = (1 - x@snps$maf)^24, get.moments = "size.based", estimation.pvalue = "kurtosis", params.sampling = list(perm.target = 100, perm.max = 50000), debug = F) }
    if(missing(H0.SKAT))  stop("'H0.SKAT should be given and can be obtained using 'NullObject.parameters()'")
    if(verbose) cat("SKAT\n")
    x.SKAT <- SKAT(x.filter, H0.SKAT, weights = SKAT.parameters$weights, maf.threshold = maf.threshold, get.moments = SKAT.parameters$get.moments, estimation.pvalue = SKAT.parameters$estimation.pvalue, params.sampling = SKAT.parameters$SKAT.sampling, debug = SKAT.parameters$debug)
  }else{
    x.SKAT <- NULL
  }
  return(list(burden = x.burden, SKAT = x.SKAT))
}