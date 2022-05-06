adjustedCADD.annotation.SNVs <- function(x, variant.scores = NULL, cores = 10, verbose = T, path.data){
  if("adjCADD" %in% colnames(x@snps)){
    warning("'adjCADD' already exists and will be replaced")
    x@snps <- x@snps[,-which(colnames(x@snps)=="adjCADD")]
  }
  if(missing(path.data)) stop("the directory 'path.data' to download and use the necessary files for RAVA-FIRST analysis should be provided")
##Check if file with scores is provided
  if(is.null(variant.scores)){
    ##Check if file with score already downloaded
    if(!file.exists(paste0(path.data, "/AdjustedCADD_v1.4_202108.tsv.gz"))){
      if(verbose){
        cat("Downloading adjusted CADD scores in ", path.data, "\n")
        curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202108.tsv.gz"), quiet = F)
      }else{
        curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202108.tsv.gz"))
      }
      curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz.tbi", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202108.tsv.gz.tbi"))
    }
    CADDfile = paste0(path.data, "/AdjustedCADD_v1.4_202108.tsv.gz")
    
    #Download README
    if(!file.exists(paste0(path.data, "/README_RAVAFIRST"))) curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/README_RAVAFIRST", destfile = paste0(path.data, "/README_RAVAFIRST"))
        
    ##Annotation with adjusted CADD
    x.pos <- paste0(x@snps$chr, ":", format(x@snps$pos-1, scientific = F, trim = T), "-", format(x@snps$pos, scientific = F, trim = T))
    ##Sort and merger positions
    x.pos.sort <- bedr.sort.region(x.pos, check.chr=F, method="natural", verbose = F)
    x.pos.merged <- bedr.merge.region(x.pos.sort, check.chr=F, verbose = F)
    ##Split pos for parallelisation
    end.indexes <- round(seq(1, length(x.pos.merged), length.out = cores+1)[-1])
    start.indexes <- c(1, end.indexes[-length(end.indexes)]+1)
    ##Get CADD score for corresponding positions
    variant.scores <- do.call(rbind, mclapply(1:cores, function(z) tabix(x.pos.merged[start.indexes[z]:end.indexes[z]], CADDfile, check.chr=F, verbose = F), mc.cores = cores))
    colnames(variant.scores) <- c("chr", "pos", "A1", "A2", "adjCADD")
  }else{
    #Annotation with provided scores file
    if(!(all(colnames(variant.scores) %in% c("chr", "pos", "A1", "A2", "adjCADD")))) stop("'variant.scores' should contain the columns 'chr', 'pos', 'A1', 'A2' and 'adjCADD'")
    if(verbose) cat("Annotation of variants with provided scores in 'variant.scores'. Warning: these scores should correspond to the adjusted CADD scores available at https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202108.tsv.gz\n")
  } 
  #Remove duplicated positions if any
  if(any(duplicated(variant.scores))){
    variant.scores <- unique(variant.scores)
  } 
  if(verbose) cat("Annotation of CADD scores \n")
  #Annotation by looking at A1->A2 and A2->A1
  x@snps$order <- 1:nrow(x@snps) #To make sure that SNPs are in the right order
  scores.notflip.tmp <- merge(x@snps, variant.scores, by = c("chr", "pos", "A1", "A2"), all.x = T)
  scores.notflip <- scores.notflip.tmp$adjCADD[order(scores.notflip.tmp$order)]
  #Inverse alleles to merge CADD scores when maf != @p
  scores.flip.tmp <- merge(x@snps, variant.scores, by.x = c("chr", "pos", "A2", "A1"), by.y = c("chr", "pos", "A1", "A2"), all.x = T)
  scores.flip <- scores.flip.tmp$adjCADD[order(scores.flip.tmp$order)]
  x@snps$adjCADD <- scores.notflip
  x@snps$adjCADD[which(is.na(x@snps$adjCADD))] <- scores.flip[which(is.na(x@snps$adjCADD))]
  x@snps <- x@snps[,-which(colnames(x@snps)=="order")]
  x
}

