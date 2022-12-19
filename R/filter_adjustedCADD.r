filter.adjustedCADD <- function(x, SNVs.scores = NULL, indels.scores = NULL, ref.level = NULL, filter = c("whole", "controls", "any"), maf.threshold = 0.01, min.nb.snps = 2, min.cumulative.maf = NULL, group = NULL, cores = 10, path.data, verbose = T){
  ##Check if bed matrix has been linked to CADD regions
  if(!("adjCADD.Median" %in% colnames(x@snps))) stop("The 'adjCADD.Median' of each CADD region should be in x@snps, please use 'set.CADDregions()'")
  #Check if x@snps$adjCADD exists and if not: annotation with CADD file
  if(!("adjCADD" %in% colnames(x@snps))){  
    ##First filter to decrease amount of data for annotation
    x <- filter.rare.variants(x, ref.level, filter, maf.threshold, min.nb.snps, min.cumulative.maf, group)
    
    x <- adjustedCADD.annotation(x, SNVs.scores = SNVs.scores, indels.scores = indels.scores, cores = cores, verbose = verbose, path.data = path.data)
  }
  else{
    if(verbose) cat("Filtering of rare variants directly on 'x@snps$adjCADD'\n")
  }
  
  #Filter CADD on median
  x <- select.snps(x, !is.na(x@snps$adjCADD) & !is.na(x@snps$adjCADD.Median) & x@snps$adjCADD >= x@snps$adjCADD.Median)
  
  #Last frequency filter
  x <- filter.rare.variants(x, ref.level, filter, maf.threshold, min.nb.snps, min.cumulative.maf, group)
  x
}
