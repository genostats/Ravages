adjustedCADD.annotation.indels <- function(x, variant.scores = NULL, cores = 10, verbose = T, path.data){
  if("adjCADD" %in% colnames(x@snps)){
    warning("'adjCADD' already exists and will be replaced")
    x@snps <- x@snps[,-which(colnames(x@snps)=="adjCADD")]
  }
  if(missing(path.data)) stop("the directory 'path.data' to download and use the necessary files for RAVA-FIRST analysis should be provided")
  if(is.null(variant.scores))  stop("Annotation of indels, 'variant.scores' should be provided with Phred scores 1.4 for Indels")
  
  if(!file.exists(paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"))){
      if(verbose){
        cat("Downloading adjusted CADD scores of indels in ", path.data, "\n")
        curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"), quiet = F)
        curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi"))
      }else{
        curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"))
        curl_download("https://lysine.univ-brest.fr/RAVA-FIRST/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi", destfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz.tbi"))
      }
    }
    
    CADDfile = paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz")
    
    ##Annotation with adjusted CADD
    x.pos <- paste0(x@snps$chr, ":", format(x@snps$pos-1, scientific = F, trim = T), "-", format(x@snps$pos, scientific = F, trim = T))
    ##Sort and merger positions
    x.pos.sort <- bedr.sort.region(x.pos, check.chr=F, verbose = F, method="natural")
    x.pos.merged <- bedr.merge.region(x.pos.sort, check.chr=F, verbose = F)
    ##Split pos for parallelisation
    end.indexes <- round(seq(1, length(x.pos.merged), length.out = cores+1)[-1])
    start.indexes <- c(1, end.indexes[-length(end.indexes)]+1)
    ##Get CADD score for corresponding positions
    tmp.scores <- do.call(rbind, mclapply(1:cores, function(z) tabix(x.pos.merged[start.indexes[z]:end.indexes[z]], CADDfile, check.chr=F, verbose = F), mc.cores = cores))
    colnames(tmp.scores) <- c("chr", "pos", "A1", "A2", "PHRED_1.4", "adjCADD", "SubRegion")
    
    x@snps$order <- 1:nrow(x@snps) #To make sure that SNPs are in the right order
    scores.notflip.tmp <- merge(x@snps, tmp.scores[,c("chr", "pos", "A1", "A2", "adjCADD")], by = c("chr", "pos", "A1", "A2"), all.x = T)
    scores.notflip <- scores.notflip.tmp$adjCADD[order(scores.notflip.tmp$order)]
    #Inverse alleles to merge CADD scores when maf != @p
    scores.flip.tmp <- merge(x@snps, tmp.scores[,c("chr", "pos", "A1", "A2", "adjCADD")], by.x = c("chr", "pos", "A2", "A1"), by.y = c("chr", "pos", "A1", "A2"), all.x = T)
    scores.flip <- scores.flip.tmp$adjCADD[order(scores.flip.tmp$order)]
    x@snps$adjCADD <- scores.notflip
    x@snps$adjCADD[which(is.na(x@snps$adjCADD))] <- scores.flip[which(is.na(x@snps$adjCADD))]
    
    ###For indels that are new: attribute the adjusted value of the nearest CADD value  
    if(sum(is.na(x@snps$adjCADD))>0){
      
      ###For indels that are new: attribute the adjusted value of the nearest CADD value  
      #Annotation of PHRED by looking at A1->A2 and A2->A1
      scores.notflip <- merge(x@snps[,c("chr", "pos", "A1", "A2", "SubRegion", "adjCADD", "order")], variant.scores, by = c("chr", "pos", "A1", "A2"), all.x = T)
      #Inverse alleles to merge CADD scores when maf != @p
      scores.flip <- merge(x@snps[,c("chr", "pos", "A1", "A2", "SubRegion", "adjCADD", "order")], variant.scores, by.x = c("chr", "pos", "A2", "A1"), by.y = c("chr", "pos", "A1", "A2"), all.x = T)
      scores.flip <- scores.flip[order(scores.flip$order),]
      tmp.x <- x@snps
      tmp.x$adjCADD <- scores.notflip$adjCADD[order(scores.notflip$order)]
      tmp.x$PHRED_1.4 <- scores.notflip$PHRED_1.4[order(scores.notflip$order)]
      tmp.x$PHRED_1.4[which(is.na(tmp.x$PHRED_1.4))] <- scores.flip[which(is.na(tmp.x$PHRED_1.4)),"PHRED_1.4"]
      tmp.x$SubRegion <- as.character(tmp.x$SubRegion)
      #Select scores of new variants
      variant.scores.new <- tmp.x[which(is.na(x@snps$adjCADD)),]
      scores.adj.indels <- fread(paste0(path.data, "/AdjustedCADD_v1.4_202204_indels.tsv.gz"), select=c(5,6,7))
      scores.byarea <- lapply(c("Coding", "Regulatory", "Intergenic"), function(z) subset(scores.adj.indels, scores.adj.indels$SubRegion==z)) 
      names(scores.byarea) <- c("Coding", "Regulatory", "Intergenic")
      #Find the nearest neighbour
      adj.scores <- mclapply(1:nrow(variant.scores.new), function(z) if(is.na(variant.scores.new[z,"PHRED_1.4"])){NA}else{as.numeric(scores.byarea[[variant.scores.new[z,"SubRegion"]]][which.min(abs(as.numeric(variant.scores.new[z,"PHRED_1.4"]) - scores.byarea[[variant.scores.new[z,"SubRegion"]]]$PHRED_1.4)), "adjCADD"])}, mc.cores = cores)
      x@snps$adjCADD[which(is.na(x@snps$adjCADD))] <- unlist(adj.scores)
      x@snps$adjCADD <- as.numeric(x@snps$adjCADD)
    }
    x@snps <- x@snps[,-which(colnames(x@snps)=="order")]
    x
}
  
    
  
