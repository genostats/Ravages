region.by.pos <- function(x, nb.groups){
  if(any(x@snps$pos == 0)) {
    stop("position equal to 0")
  }
  x <- select.snps(x, is.autosome(x))
  dist <- c(abs(diff(x@snps$pos)),0)
  pos_d <- data.frame(chr=x@snps$chr,id=x@snps$id,pos=x@snps$pos,d=dist)
       
  seuil <- 1000
  while(length(which(pos_d$d>seuil)) > nb.groups){
    seuil <- seuil+1000  }
                     
  diff_pos <- diff(c(0,which(pos_d$d>seuil), length(pos_d$d)))
  gpe <- as.factor(rep.int(1:length(diff_pos), diff_pos))  
                           
  gpe_geno <- data.frame(pos_d,gpe=gpe)

  for (i in 1:21){
  if(tail(gpe_geno$gpe[which(gpe_geno$chr==i)],1)==head(gpe_geno$gpe[which(gpe_geno$chr==i+1)],1)){
    warnings("SNPs from 2 different chromosomes in the same group")}} 

  return(gpe_geno)
}

