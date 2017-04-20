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
                     
  diff_pos <- diff(c(0,which(pos_d$d>seuil), length(pos_d$d)))
  gpe <- as.factor(rep.int(1:length(diff_pos), diff_pos))  
                           
  gpe_geno <- data.frame(pos_d,gpe=gpe)

  for (i in 1:21){
  if(tail(gpe_geno$gpe[which(gpe_geno$chr==i)],1)==head(gpe_geno$gpe[which(gpe_geno$chr==i+1)],1)){
    warnings("SNPs from 2 different chromosomes in the same group")}} 

  return(gpe_geno)
}



B_cast <- function(x, regions, threshold=0.01){

  weights.0 <- ifelse( (x@snps$maf == (1-x@p) & x@snps$maf < threshold), 1, 0)
  weights.1 <- ifelse( (x@snps$maf < threshold), 1, 0)
  weights.2 <- ifelse( (x@snps$maf == x@p & x@snps$maf < threshold), 1, 0)
      
  if(length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(regions) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  if(!is.factor(regions)) stop("'regions' is not a factor")
                
  B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(regions), regions, weights.0, weights.1, weights.2)
  colnames(B) <- levels(regions)
  rownames(B) <- x@ped$id
  Bcast <- (B>0)+0
  return(Bcast)
  }


Chi2 <- function(centre, B){
 x <- table(centre,B)
 if (is.matrix(x)) {
   nr <- as.integer(nrow(x))
   nc <- as.integer(ncol(x))
   n <- sum(x)
     if (is.na(nr) || is.na(nc) || is.na(nr * nc))
       stop("invalid nrow(x) or ncol(x)", domain = NA)
   sr <- rowSums(x)
   sc <- colSums(x)
   E <- outer(sr, sc, "*")/n
                                    
   if (all(E>5)) {
     c <- chisq.test( x )
   }
   else {
     c <- chisq.test( x, simulate.p.value=TRUE, B=1e5)
   }
 }
 return(c$p.value)       
 }



WSS <- function(x, regions){

 Q <- (2*x@snps$N2+x@snps$N1 + 1) / ( 2*(x@snps$N0+x@snps$N1+x@snps$N2) +2 )
 W <- sqrt((nrow(x)-x@snps$NAs) * Q * (1-Q))
 weights.0 <- ifelse( Q > 0.5, 2/W, 0)
 weights.1 <- 1/W
 weights.2 <- ifelse( Q < 0.5, 2/W, 0)
    
 if(length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(regions) != ncol(x)) {
   stop("x and weights dimensions mismatch")
   }
 if(!is.factor(regions)) stop("'regions' is not a factor")
                            
   B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(regions), regions, weights.0, weights.1, weights.2)
   colnames(B) <- levels(regions)
   rownames(B) <- x@ped$id
   return(B)
 }



burden <- function(x, regions, weights.0, weights.1, weights.2) {
  if(length(weights.0) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(regions) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  if(!is.factor(regions)) stop("'regions' is not a factor")

  B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(regions), regions, weights.0, weights.1, weights.2)
  colnames(B) <- levels(regions)
  rownames(B) <- x@ped$id
  B
}
