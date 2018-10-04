region.by.gene.old <- function(x, genes, include.all=FALSE){
  #Remove gene on different chromosomes
  if(any(table(genes$Gene_Name)!=1)==TRUE){
    a <- unique(as.character(genes$Gene_Name))
    genes$Gene_Name <- factor(as.character(genes$Gene_Name), levels=a)
    w <- names(which(table(genes$Gene_Name)==1))
    genes <- subset(genes, genes$Gene_Name %in% w)
    genes$Gene_Name <- droplevels(genes$Gene_Name)
  } 
  
  #Remove gap between the genes
  if(include.all==TRUE){
    genes_c <- genes
    for (i in 2:nrow(genes)){
      if(genes$Start[i] > genes$End[i-1]){
        genes_c$Start[i] <- mean(c(genes$End[i-1], genes$Start[i]))
        genes_c$End[i-1] <- genes_c$Start[i]
      }
    }  
  genes <- genes_c  
  }
  
  gpe <- data.frame(chr=x@snps$chr, pos=x@snps$pos, id=x@snps$id, gpe=factor( rep(NA, nrow(x@snps)), levels = levels(genes$Gene_Name)) )
  for (i in 1:nrow(genes)){
    w <- (x@snps$chr == genes$Chr[i] & x@snps$pos >= genes$Start[i] & x@snps$pos <= genes$End[i])
    gpe$gpe[w] <- genes$Gene_Name[i]
  }
  gpe$gpe <- droplevels(gpe$gpe)
  return(gpe)
  }
  


region.by.gene <- function(x, genes, include.all = FALSE) {
  # remove duplicated genes if any
  w <- duplicated(genes$Gene_Name)
  if(any(w)) {
    genes <- genes[!w,]
  }

  # check if genes is sorted by chr / starting pos
  n <- nrow(genes)
  chr1 <- genes$Chr[1:(n-1)]
  chr2 <- genes$Chr[2:n]
  b <- (chr1 < chr2) | (chr1 == chr2 & genes$Start[1:(n-1)] <= genes$Start[2:n])
  if(!all(b)) {
    genes <- genes[ order(genes$Chr, genes$Start), ]
  }
  
  # if include.all, define larger regions
  if(include.all) {
    M <- as.integer(max(x@snps$pos, genes$End))  # joue le rôle de la position infinie !
    b <- genes$Chr[2:n] == genes$Chr[1:(n-1)]
    start <- ifelse(b, as.integer(0.5*(genes$Start[-1] + genes$End[-1])),  0L)
    end <- ifelse(b, start-1L, M)
    genes$Start <- c(0L,start)
    genes$End <- c(end,M)
  }
  
  gpe <- data.frame(chr = x@snps$chr, pos = x@snps$pos, id = x@snps$id, gpe = factor( rep(NA, nrow(x@snps)), levels = levels(genes$Gene_Name)) )
 
  R <- .Call("oz_label_genes", PACKAGE = "Ravages", genes$Chr, genes$Start, genes$End, x@snps$chr, x@snps$pos)
  R <- ifelse(R == 0, NA, R)
  levels(R) <- as.character(genes$Gene_Name)
  class(R) <- "factor"

  # remettre les niveaux dans l'ordre de ceux de genes$Gene_Name
  R <- factor(R, levels(genes$Gene_Name))

  gpe$gpe <- R

  gpe
}


