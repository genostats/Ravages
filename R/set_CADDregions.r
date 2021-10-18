set.CADDregions <- function(x, verbose = T) {
  Ravages_path <- system.file(package="Ravages")
  ##Check if file with score already downloaded
  if(!file.exists(paste0(Ravages_path, "/README_RAVAFIRST"))) curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/README_RAVAFIRST", destfile = paste0(Ravages_path, "/README_RAVAFIRST")) #download README
  if(!file.exists(paste0(Ravages_path, "/CADDRegions.2021.hg19.bed.gz"))){
    if(verbose) cat("Downloading CADD regions in ", Ravages_path, "\n")
    curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/CADDRegions.2021.hg19.bed.gz", destfile = paste0(Ravages_path, "/CADDRegions.2021.hg19.bed.gz"))
    curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/CADDRegions.2021.hg19.bed.gz.tbi", destfile = paste0(Ravages_path, "/CADDRegions.2021.hg19.bed.gz.tbi"))
  }
  if(!file.exists(paste0(Ravages_path, "/FunctionalAreas.hg19.bed.gz"))){
    if(verbose) cat("Downloading Genomic Categories in ", Ravages_path, "\n")
    curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/FunctionalAreas.hg19.bed.gz", destfile = paste0(Ravages_path, "/FunctionalAreas.hg19.bed.gz"))
    curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/FunctionalAreas.hg19.bed.gz.tbi", destfile = paste0(Ravages_path, "/FunctionalAreas.hg19.bed.gz.tbi"))
  }

  regions <- read.table(gzfile(paste0(Ravages_path, "/CADDRegions.2021.hg19.bed.gz")), header = T, as.is = T)
  regions$Name <- factor(regions$Name, levels = unique(regions$Name))
  subregions <- read.table(gzfile(paste0(Ravages_path, "/FunctionalAreas.hg19.bed.gz")), header = T, as.is = T)

  #Add one to start to take into account bed format
  regions$Start <- regions$Start+1
  subregions$Start <- subregions$Start+1

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






