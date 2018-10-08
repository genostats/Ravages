set.genomic.region.by.gene <- function(x, genes, include.all = FALSE) {
  # ce test est OK pour les facteurs aussi
  if(typeof(x@snps$chr) != "integer") 
    stop("x@snps$chr should be either a vector of integers, or a factor with same levels as genes$Chr")
  
  # difficile de tester intelligemment l'identité des niveaux ...
  if(is.factor(x@snps$chr)) 
    genes$Chr <- as.factor(genes$Chr)
   
  
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
  
  R <- .Call("oz_label_genes", PACKAGE = "Ravages", genes$Chr, genes$Start, genes$End, x@snps$chr, x@snps$pos)
  R <- ifelse(R == 0, NA, R)
  levels(R) <- as.character(genes$Gene_Name)
  class(R) <- "factor"

  # remettre les niveaux dans l'ordre de ceux de genes$Gene_Name 
  R <- factor(R, levels(genes$Gene_Name))

  x@snps$genomic.region <- droplevels(R) ## attention au droplevels !

  x
}


