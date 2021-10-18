RAVA.FIRST <- function(x, variant.scores = NULL, ref.level, filter = c("whole", "controls", "any"), maf.threshold = 0.01, min.nb.snps = 2, min.cumulative.maf = NULL, group = NULL, cores = 10, burden = TRUE, H0.burden, burden.parameters, SKAT = TRUE, H0.SKAT, SKAT.parameters, verbose = TRUE){
  ##Check if good paramters for RVAT are given
  if(burden & missing(H0.burden)) stop("'H0.burden should be given if 'burden = TRUE' and can be obtained using 'NullObject.parameters()'")
  if(SKAT & missing(H0.SKAT))  stop("'H0.SKAT should be given if 'SKAT = TRUE' and can be obtained using 'NullObject.parameters()'")
  
  ###Attribute CADD regions to variants
  if(verbose) cat("Attribution of CADD regions\n")
  x <- set.CADDregions(x)
  
  #Importation of CADD regions to give position information in final results
  Ravages_path <- system.file(package="Ravages")
  regions <- read.table(gzfile(paste0(Ravages_path, "/CADDRegions.2021.hg19.bed.gz")), header = T, as.is = T)
  rownames(regions) <- regions$Name
  
  ###Filtering based on CADD scores
  x.filter <- filter.adjustedCADD(x, variant.scores = variant.scores, ref.level = ref.level, filter = filter, maf.threshold = maf.threshold, min.nb.snps = min.nb.snps, min.cumulative.maf = min.cumulative.maf, group = group, verbose = verbose)
  if(verbose) cat("Filtering of rare variants within CADD regions\n")
  ###Association
  if(burden){
    if(missing(burden.parameters)){ burden.parameters <- list(burden.function = WSS, get.effect.size = F) }
    if(verbose) cat("Burden test\n")
    x.burden <- burden.subscores(x.filter, H0.burden, burden.function = burden.parameters$burden.function, maf.threshold = maf.threshold, get.effect.size = burden.parameters$get.effect.size, alpha = 0.05, cores = cores, verbose = verbose)
    #Add information about CADD region
    x.burden <- cbind(x.burden, regions[rownames(x.burden),c(1:3,5,7)])
  }else{
    x.burden <- NULL
  }
  if(SKAT){
    if(missing(SKAT.parameters)){ SKAT.parameters <- list(get.moments = "size.based", estimation.pvalue = "kurtosis", params.sampling = list(perm.target = 100, perm.max = 50000), debug = F) }
    if(verbose) cat("SKAT\n")
    x.SKAT <- SKAT(x.filter, H0.SKAT, weights = (1 - x.filter@snps$maf)^24, maf.threshold = maf.threshold, get.moments = SKAT.parameters$get.moments, estimation.pvalue = SKAT.parameters$estimation.pvalue, params.sampling = SKAT.parameters$params.sampling, debug = SKAT.parameters$debug)
    #Add information about CADD region
    x.SKAT <- cbind(x.SKAT, regions[rownames(x.SKAT),c(1:3,5,7)])
  }else{
    x.SKAT <- NULL
  }
  
  return(list(burden = x.burden, SKAT = x.SKAT))
}
