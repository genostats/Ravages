Beta.M.2 <- function(x, centre, region, which.snps = rep(TRUE, ncol(x)), target = 0, B.max = 0) {
  r <- .Call('oz_beta_m', PACKAGE = 'oz', x@bed, which.snps, region, centre, target, B.max) 
  res <- data.frame( genomic.region = levels(region), r )
  if(!is.null(res$A)) res$p <- (res$A+1)/(res$B+1)
  return(res)
}
