adjustedCADD.annotation <- function(x, cores = 10, verbose = T){
  Ravages_path <- system.file(package="Ravages")

  ##Check if file with score already downloaded
  if(!file.exists(paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz"))){
    if(verbose){
      cat("Downloading adjusted CADD scores in ", Ravages_path, "\n")
      curl_download("https://lysine.univ-brest.fr/CADD_Regions/AdjustedCADD_v1.4_202108.tsv.gz", destfile = paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz"), quiet = F)
    }else{
      curl_download("https://lysine.univ-brest.fr/CADD_Regions/AdjustedCADD_v1.4_202108.tsv.gz", destfile = paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz"))
    }
    curl_download("https://lysine.univ-brest.fr/CADD_Regions/AdjustedCADD_v1.4_202108.tsv.gz.tbi", destfile = paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz.tbi"))
  }
  CADDfile = paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz")
    
  ##Annotation with adjusted CADD
  x.pos <- paste0(x@snps$chr, ":", format(x@snps$pos-1, scientific = F, trim = T), "-", format(x@snps$pos, scientific = F, trim = T))
  ##Sort and merger positions
  x.pos.sort <- bedr.sort.region(x.pos, check.chr=F, verbose = F)
  x.pos.merged <- bedr.merge.region(x.pos.sort, check.chr=F, verbose = F)
  ##Split pos for parallelisation
  end.indexes <- round(seq(1, length(x.pos.merged), length.out = cores+1)[-1])
  start.indexes <- c(1, end.indexes[-length(end.indexes)]+1)
  ##Get CADD score for corresponding positions
  if(verbose) cat("Annotation of CADD scores \n")
  x.scores <- do.call(rbind, mclapply(1:cores, function(z) tabix(x.pos.merged[start.indexes[z]:end.indexes[z]], CADDfile, check.chr=F, verbose = F), mc.cores = cores))
  colnames(x.scores) <- c("chr", "pos", "A1", "A2", "adjCADD")
  #Annotation by looking at A1->A2 and A2->A1
  scores.notflip <- merge(x@snps, x.scores, by = c("chr", "pos", "A1", "A2"), sort.x = F, all.x = T)$adjCADD
  #Inverse alleles to merge CADD scores when maf != @p
  scores.flip <- merge(x@snps, x.scores, by.x = c("chr", "pos", "A2", "A1"), by.y = c("chr", "pos", "A1", "A2"), sort.x = F, all.x = T)$adjCADD
  x@snps$adjCADD <- scores.notflip
  x@snps$adjCADD[which(is.na(x@snps$adjCADD))] <- scores.flip[which(is.na(x@snps$adjCADD))]
  x
}
