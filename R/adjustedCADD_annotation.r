adjustedCADD.annotation <- function(x, SNVs.scores = NULL, indels.scores = NULL, cores = 10, verbose = T, path.data){
  if(missing(path.data)) stop("the directory 'path.data' to download and use the necessary files for RAVA-FIRST analysis should be provided")
  if(any(nchar(x@snps$A1)>1) | any(nchar(x@snps$A2)>1)){
    if(is.null(indels.scores))  stop("Indels are present, 'indels.scores' should be provided with Phred scores 1.4 for Indels")
    if(!("SubRegion" %in% colnames(x@snps))) stop("The 'SubRegion' of each CADD region should be in x@snps, please use 'set.CADDregions()'")
    which.indels <- (nchar(x@snps$A1)>1 | nchar(x@snps$A2)>1)
    x.indels <- select.snps(x, which.indels)
    if(ncol(x.indels) != nrow(indels.scores)) stop("'indels.scores' has wrong dimensions")
    x.indels <- adjustedCADD.annotation.indels(x.indels, variant.scores = indels.scores, cores = cores, verbose = verbose, path.data = path.data)
    if(ncol(x.indels)<ncol(x)){
      x.SNVs <- select.snps(x, !which.indels)
      x.SNVs <- adjustedCADD.annotation.SNVs(x.SNVs, variant.scores = SNVs.scores, cores = cores, verbose = verbose, path.data = path.data)
      x@snps$adjCADD[which.indels] <- x.indels@snps$adjCADD
      x@snps$adjCADD[!which.indels] <- x.SNVs@snps$adjCADD
    }else{
      x <- x.indels
    }
  }else{
    x <- adjustedCADD.annotation.SNVs(x, variant.scores = SNVs.scores, cores = cores, verbose = verbose, path.data = path.data)
  }
  x
}



