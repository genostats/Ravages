set.CADDregions <- function(x, verbose = T) {
  Ravages_path <- system.file(package="Ravages")
  ##Check if file with score already downloaded
  if(!file.exists(paste0(Ravages_path, "/CADDRegions.2021.hg19.tsv.gz"))){
    if(verbose) cat("Downloading CADD regions in ", Ravages_path, "\n")
    curl_download("https://lysine.univ-brest.fr/CADD_Regions/CADDRegions.2021.hg19.tsv.gz", destfile = paste0(Ravages_path, "/CADDRegions.2021.hg19.tsv.gz"))
  }
  if(!file.exists(paste0(Ravages_path, "/FunctionalAreas.hg19.tsv.gz"))){
    if(verbose) cat("Downloading Genomic Areas in ", Ravages_path, "\n")
    curl_download("https://lysine.univ-brest.fr/CADD_Regions/FunctionalAreas.hg19.tsv.gz", destfile = paste0(Ravages_path, "/FunctionalAreas.hg19.tsv.gz"))
  }

  regions <- read.table(gzfile(paste0(Ravages_path, "/CADDRegions.2021.hg19.tsv.gz")), header = T, as.is = T)
  regions$Name <- factor(regions$Name, levels = unique(regions$Name))
  subregions <- read.table(gzfile(paste0(Ravages_path, "/FunctionalAreas.hg19.tsv.gz")), header = T, as.is = T)

  R <- .Call("label_multiple_genes", PACKAGE = "Ravages", regions$Chr, regions$Start, regions$End, x@snps$chr, x@snps$pos)
  R.genename <- unlist(lapply(R, function(z) paste(levels(regions$Name)[unlist(z)], collapse=",")))  
  R.genename[which(R.genename=="")] <- NA
  #############Add median CADD observed by genomic region
  R.median <- unlist(lapply(R, function(z) paste(regions$Median[unlist(z)], collapse=",")))  
  R.median[which(R.median=="")] <- NA

  x@snps$genomic.region <- R.genename
  x@snps$genomic.region <- factor(x@snps$genomic.region, levels = unique(x@snps$genomic.region))
  x@snps$adjCADD.Median <- R.median
  
  ##################Add subregions
  Rsub <- .Call("label_multiple_genes", PACKAGE = "Ravages", subregions$Chr, subregions$Start, subregions$End, x@snps$chr, x@snps$pos)
  Rsub.genename <- unlist(lapply(Rsub, function(z) paste(subregions$GenomicArea[unlist(z)], collapse=",")))  
  Rsub.genename[which(Rsub.genename=="")] <- NA

  x@snps$SubRegion <- Rsub.genename 
  
  x
}






