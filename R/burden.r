region_geno <- function(x, nb_groupes){
  if(any(x@snps$pos == 0)) {
    stop("position equal to 0")
  }
  
  dist <- c(abs(diff(x@snps$pos)),0)
  pos_d <- data.frame(chr=x@snps$chr,id=x@snps$id,pos=x@snps$pos,d=dist)
  pos_d <- pos_d[which(pos_d$chr!=23 & pos_d$chr!=24 & pos_d$chr!=26),]
       
  seuil <- 1000
  while(length(which(pos_d$d>seuil)) > nb_groupes){
    seuil <- seuil+1000  }
                     
  diff_pos <- diff(c(0,which(pos_d$chr!=23 & pos_d$chr!=24 & pos_d$chr!=26 & pos_d$d>seuil)))
  gpe <- as.factor(rep.int(1:length(diff_pos), diff_pos))  
                           
  gpe_geno <- data.frame(pos_d,gpe=gpe)

  for (i in 1:21){
  if(tail(gpe_geno$gpe[which(gpe_geno$chr==i)],1)==head(gpe_geno$gpe[which(gpe_geno$chr==i+1)],1)){
    warnings("SNPs from 2 different chromosomes in the same group")}} 

  return(gpe_geno)
}



B_cast <- function(x, seuil, regions){

  weights.0 <- ifelse( (x@snps$maf == (1-x@p) & x@snps$maf<seuil), 1, 0)
  weights.1 <- ifelse( (x@snps$maf<0.01), 1, 0)
  weights.2 <- ifelse( (x@snps$maf == x@p & x@snps$maf<0.01), 1, 0)
      
  if(length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(regions) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  if(!is.factor(regions)) stop("'regions' is not a factor")
                
  B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(regions), regions, weights.0, weights.1, weights.2)
  colnames(B) <- levels(regions)
  rownames(B) <- x@ped$id
  Bcast <- (B>0)+0
  return(Bcast)
  }



burden <- function(x, regions, weights.0, weights.1, weights.2) {
  if(length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(regions) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  if(!is.factor(regions)) stop("'regions' is not a factor")

  B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(regions), regions, weights.0, weights.1, weights.2)
  colnames(B) <- levels(regions)
  rownames(B) <- x@ped$id
  B
}

