set.genomic.region.subregion <- function(x, regions, subregions, split = TRUE) {

  if(!is.factor(regions$Name) | !is.factor(subregions$Name)) stop("regions$Name and subregions$Name should be factors with levels ordered in the genome order")

  # ce test est OK pour les facteurs aussi
  if(typeof(x@snps$chr) != "integer") 
    stop("x@snps$chr should be either a vector of integers, or a factor with same levels as regions$Chr")
  
  #Add one to start to take into account bed format
  regions$Start <- regions$Start + 1 
  
  # remove duplicated regions if any
  w <- duplicated(regions$Name)
  if(any(w)) {
    regions <- regions[!w,]
  }

  # check if regions is sorted by chr / starting pos
  n <- nrow(regions)
  chr1 <- regions$Chr[1:(n-1)]
  chr2 <- regions$Chr[2:n]
  b <- (chr1 < chr2) | (chr1 == chr2 & regions$Start[1:(n-1)] <= regions$Start[2:n])
  if(!all(b)) {
    regions <- regions[ order(regions$Chr, regions$Start), ]
  }
  
  # check if subregions is sorted by chr / starting pos
  n <- nrow(subregions)
  chr1 <- subregions$Chr[1:(n-1)]
  chr2 <- subregions$Chr[2:n]
  b <- (chr1 < chr2) | (chr1 == chr2 & subregions$Start[1:(n-1)] <= subregions$Start[2:n])
  if(!all(b)) {
    subregions <- subregions[ order(subregions$Chr, subregions$Start), ]
  }
  
  ###Annotation region
  R <- .Call("label_multiple_genes", PACKAGE = "Ravages", regions$Chr, regions$Start, regions$End, x@snps$chr, x@snps$pos)
  R.genename <- unlist(lapply(R, function(z) paste(levels(regions$Name)[unlist(z)], collapse=",")))  
  R.genename[which(R.genename=="")] <- NA

  x@snps$genomic.region <- R.genename
  x@snps$genomic.region <- factor(x@snps$genomic.region, levels = unique(x@snps$genomic.region))
  if(any(grepl(x@snps$genomic.region, pattern = ","))){
    x <- bed.matrix.split.genomic.region(x, genomic.region = x@snps$genomic.region, split.pattern = ",")
  }
  
  ###Annotation subregion
  Rsub <- .Call("label_multiple_genes", PACKAGE = "Ravages", subregions$Chr, subregions$Start, subregions$End, x@snps$chr, x@snps$pos)
  Rsub.genename <- unlist(lapply(Rsub, function(z) paste(subregions$Name[unlist(z)], collapse=",")))  
  Rsub.genename[which(Rsub.genename=="")] <- NA

  x@snps$SubRegion <- Rsub.genename 
  x@snps$SubRegion <- factor(x@snps$SubRegion, levels = unique(x@snps$SubRegion))
  #if(any(grepl(x@snps$SubRegion, pattern = ","))){
  #  x <- bed.matrix.split.genomic.region(x, genomic.region = x@snps$SubRegion, split.pattern = ",")
  #}
  
  x
}


