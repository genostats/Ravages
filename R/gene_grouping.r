pos_gpe <- function(x, nb_groups){
  if(any(x@snps$pos == 0)) {
    stop("position equal to 0")
  }
  
  dist <- c(abs(diff(x@snps$pos)),0)
  pos_d <- data.frame(chr=x@snps$chr,id=x@snps$id,pos=x@snps$pos,d=dist)
  pos_d <- pos_d[which(pos_d$chr!=23 & pos_d$chr!=24 & pos_d$chr!=26),]
       
  seuil <- 1000
  while(length(which(pos_d$d>seuil)) > nb_groups){
    seuil <- seuil+1000  }
                     
  diff_pos <- diff(c(0,which(pos_d$d>seuil), length(pos_d$d)))
  gpe <- as.factor(rep.int(1:length(diff_pos), diff_pos))  
                           
  gpe_geno <- data.frame(pos_d,gpe=gpe)

  for (i in 1:21){
  if(tail(gpe_geno$gpe[which(gpe_geno$chr==i)],1)==head(gpe_geno$gpe[which(gpe_geno$chr==i+1)],1)){
    warnings("SNPs from 2 different chromosomes in the same group")}} 

  return(gpe_geno)
}



genes_gpe <- function(x, genes, include.all=FALSE){
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
  

