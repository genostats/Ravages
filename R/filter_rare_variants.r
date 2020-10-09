
filter.rare.variants <- function(x, ref.level, filter = c("whole", "controls", "any"), maf.threshold = 0.01, min.nb.snps, group) {
  if(!is.factor(x@snps$genomic.region)) stop("'x@snps$genomic.region' should be a factor")
  x@snps$genomic.region <- droplevels(x@snps$genomic.region)

  #Check if good filter
  filter <- match.arg(filter)

  if(filter != "whole"){  
    if(missing(group)){
      group <- x@ped$pheno
    }else{
      #Check dim of group
      if(length(group) != nrow(x@ped)) stop("group has wrong length")
    }
    if(!is.factor(group)) stop("'group' should be a factor")
  }

  filter <- match.arg(filter)
  if(filter == "controls") {
    if(missing(ref.level)) stop("Need to specify the controls group")
    which.controls <- group == ref.level
    st <- .Call('gg_geno_stats_snps', PACKAGE = "Ravages", x@bed, rep(TRUE, ncol(x)), which.controls)$snps
    p <- (st$N0.f + 0.5*st$N1.f)/(st$N0.f + st$N1.f + st$N2.f)
    maf <- pmin(p, 1-p)
    w <- (maf < maf.threshold)
  }
  if(filter == "any"){
    # filter = any
    w <- rep(FALSE, ncol(x))
    for(i in unique(group)) {
      which.c <- (group == i)
      st <- .Call('gg_geno_stats_snps', PACKAGE = "Ravages", x@bed, rep(TRUE, ncol(x)), which.c)$snps
      p <- (st$N0.f + 0.5*st$N1.f)/(st$N0.f + st$N1.f + st$N2.f)
      maf <- pmin(p, 1-p)
      w <- w | (maf < maf.threshold)
    }
  }
  if(filter == "whole"){
      ##Filter = whole
      w <- (x@snps$maf < maf.threshold)
  }

  x <- select.snps(x, w & x@snps$maf > 0)
  if(!missing(min.nb.snps)) {
    nb.snps <- table(x@snps$genomic.region)
    keep <- names(nb.snps)[nb.snps >= min.nb.snps]
    x <- select.snps(x, x@snps$genomic.region %in% keep)
  }
    x@snps$genomic.region <- droplevels(x@snps$genomic.region)
  x
}
