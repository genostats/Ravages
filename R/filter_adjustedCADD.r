filter.adjustedCADD <- function(x, variant.scores = NULL, ref.level, filter = c("whole", "controls", "any"), maf.threshold = 0.01, min.nb.snps = 2, min.cumulative.maf, group, cores = 10, verbose = T) {
  ##Check if bed matrix has been linked to CADD regions
    if(!("adjCADD.Median" %in% colnames(x@snps))) stop("The 'adjCADD.Median' of each CADD region should be in x@snps, please use 'set.CADDregions()'")
  #Check if x@snps$adjCADD exists
  if(!("adjCADD" %in% colnames(x@snps))){  
    Ravages_path <- system.file(package="Ravages")
    
  
    ##First filter to decrease amount of data for annotation
    x <- filter.rare.variants(x, ref.level, filter, maf.threshold, min.nb.snps, min.cumulative.maf, group)
  
    ##Check if file with scores is provided
    if(is.null(variant.scores)){
      ##Check if file with score already downloaded
      if(!file.exists(paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz"))){
        if(verbose){
          cat("Downloading adjusted CADD scores in ", Ravages_path, "\n")
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz", destfile = paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz"), quiet = F)
        }else{
          curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz", destfile = paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz"))
        }
        curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz.tbi", destfile = paste0(Ravages_path, "/AdjustedCADD_v1.4_202108.tsv.gz.tbi"))
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
      variant.scores <- do.call(rbind, mclapply(1:cores, function(z) tabix(x.pos.merged[start.indexes[z]:end.indexes[z]], CADDfile, check.chr=F, verbose = F), mc.cores = cores))
      colnames(variant.scores) <- c("chr", "pos", "A1", "A2", "adjCADD")
    }else{
      #Annotation with provided scores file
      if(!(all(colnames(variant.scores) %in% c("chr", "pos", "A1", "A2", "adjCADD")))) stop("'variant.scores' should contain the columns 'chr', 'pos', 'A1', 'A2' and 'adjCADD'")
      if(verbose) cat("Annotation of variants with provided scores in 'variant.scores'. Warning: these scores should correspond to the adjusted CADD scores available at https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz\n")
    }
    ##Get CADD score for corresponding positions
    if(verbose) cat("Annotation of CADD scores \n")
    #Annotation by looking at A1->A2 and A2->A1
    scores.notflip <- merge(x@snps, variant.scores, by = c("chr", "pos", "A1", "A2"), sort.x = F, all.x = T)$adjCADD
    #Inverse alleles to merge CADD scores when maf != @p
    scores.flip <- merge(x@snps, variant.scores, by.x = c("chr", "pos", "A2", "A1"), by.y = c("chr", "pos", "A1", "A2"), sort.x = F, all.x = T)$adjCADD
    x@snps$adjCADD <- scores.notflip
    x@snps$adjCADD[which(is.na(x@snps$adjCADD))] <- scores.flip[which(is.na(x@snps$adjCADD))]
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
