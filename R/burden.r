B_cast <- function(x, genomic_region, maf.threshold=0.01){

  weights.0 <- ifelse( (x@snps$maf == (1-x@p) & x@snps$maf < maf.threshold), 1, 0)
  weights.1 <- ifelse( (x@snps$maf < maf.threshold), 1, 0)
  weights.2 <- ifelse( (x@snps$maf == x@p & x@snps$maf < maf.threshold), 1, 0)
      
  if(length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(genomic_region) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  if(!is.factor(genomic_region)) stop("'genomic_region' is not a factor")
                
  B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(genomic_region), genomic_region, weights.0, weights.1, weights.2)
  colnames(B) <- levels(genomic_region)
  rownames(B) <- x@ped$id
  Bcast <- (B>0)+0
  return(Bcast)
  }


Chi2 <- function(B_CAST, group){
 x <- table(group, B_CAST)
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



WSS <- function(x, genomic_region){

 Q <- (2*x@snps$N2+x@snps$N1 + 1) / ( 2*(x@snps$N0+x@snps$N1+x@snps$N2) +2 )
 W <- sqrt((nrow(x)-x@snps$NAs) * Q * (1-Q))
 weights.0 <- ifelse( Q > 0.5, 2/W, 0)
 weights.1 <- 1/W
 weights.2 <- ifelse( Q < 0.5, 2/W, 0)
    
 if(length(weights.0) != ncol(x) | length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(genomic_region) != ncol(x)) {
   stop("x and weights dimensions mismatch")
   }
 if(!is.factor(genomic_region)) stop("'genomic_region' is not a factor")
                            
   B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(genomic_region), genomic_region, weights.0, weights.1, weights.2)
   colnames(B) <- levels(genomic_region)
   rownames(B) <- x@ped$id
   return(B)
 }

