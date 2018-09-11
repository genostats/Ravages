region.by.gene <- function(x, genes, include.all=FALSE){
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
  

