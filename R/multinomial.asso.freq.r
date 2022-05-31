multinomial.asso.freq <- function(x, pheno = x@ped$pheno, ref.level, test = c("Genotypic", "Allelic"), get.effect.size = F, min.maf.threshold = 0.05){
  test <- match.arg(test)
  if(!(test %in% c("Genotypic", "Allelic"))) stop("Only 'Genotypic' or 'Allelic' chi-square 'test' are available")
  if (!is.factor(pheno))  stop("'pheno' is not a factor")
  pheno <- droplevels(pheno)
  if(get.effect.size)
    if(is.null(ref.level)) stop("'ref.level' should be specified to estimate OR")
  #Keep only frequent variants
  x <- select.snps(x, x@snps$maf >= min.maf.threshold)
  if(get.effect.size) pheno <- relevel(pheno, ref = ref.level)
  x@ped$pheno <- pheno
  if(test == "Genotypic"){
    counts.bygroups <- do.call(cbind, sapply(levels(pheno), function(z) select.inds(x, pheno == z)@snps[,c("N0", "N1", "N2")]))
    pval <- data.frame(chr = x@snps$chr, pos = x@snps$pos, p.value = sapply(1:nrow(counts.bygroups), function(z) suppressWarnings(chisq.test(matrix(as.numeric(counts.bygroups[z,]), nrow = nlevels(pheno), ncol = 3, byrow = T))$p.value)))
    rownames(pval) <- x@snps$id
    #Compute OR: one for each genotype comparing to A1A1
    if(get.effect.size){ 
      OR <- lapply( 1:(nlevels(pheno)-1), function(gpe){ res.OR <- t(sapply(1:nrow(counts.bygroups), function(z) matrix(c((counts.bygroups[z,1]*counts.bygroups[z,3*gpe+2])/(counts.bygroups[z,2]*counts.bygroups[z,3*gpe+1]), (counts.bygroups[z,1]*counts.bygroups[z,3*gpe+3])/(counts.bygroups[z,3]*counts.bygroups[z,3*gpe+1])), ncol = 2))) ; rownames(res.OR) <- x@snps$id ; colnames(res.OR) <- c("Hetero", "Homo_alt") ; return(res.OR) } )
      names(OR) <- levels(pheno)[-1]
    }
  }
  if(test == "Allelic"){
    counts.bygroups <- do.call(cbind, sapply(levels(pheno), function(z){ tmp <- select.inds(x, pheno == z)@snps[,c("N0", "N1", "N2")] ; c(tmp[1]*2+tmp[2], tmp[2]+tmp[3]*2)}))
    pval <- data.frame(chr = x@snps$chr, pos = x@snps$pos, p.value = sapply(1:nrow(counts.bygroups), function(z) suppressWarnings(chisq.test(matrix(as.numeric(counts.bygroups[z,]), nrow = nlevels(pheno), ncol = 2, byrow = T))$p.value)))
    rownames(pval) <- x@snps$id
    #Compute OR: allele A2 vs A1
    if(get.effect.size){
      OR <- lapply( 1:(nlevels(pheno)-1), function(gpe){ res.OR <- sapply(1:nrow(counts.bygroups), function(z) (counts.bygroups[z,1]*counts.bygroups[z,2*gpe+2])/(counts.bygroups[z,2]*counts.bygroups[z,2*gpe+1])) ; names(res.OR) <- x@snps$id ; return(res.OR) })
      names(OR) <- levels(pheno)[-1]
    }
  }
  if(get.effect.size) res <- list(Asso  = pval, OR = OR)
  else res <- pval
  res
}
  
