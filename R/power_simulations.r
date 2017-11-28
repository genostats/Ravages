

# pas inspiré pour trouver un nom
fff <- function(pop.maf, size, baseline, replicates, OR.pars) {
  OR.pars$n.variants <- length(pop.maf)
  ff <- function() {
    OR <- do.call( OR.matrix, OR.pars)
    MAFS <- group.mafs(pop.maf, OR, baseline)
    random.bed.matrix(MAFS, size)
  }
  x <- suppressWarnings(do.call( cbind, replicate(replicates, ff())))  # warnings à cause 
  x <- set.stats(x, verbose = FALSE)
  x@snps$genomic.region <- factor( rep(sprintf("R%d", 1:replicates), each = length(pop.maf)) )
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x
}

filter.rare.variants <- function(x, filter = c("controls", "any"), maf.threshold = 0.01, min.nb.snps) {
  filter <- match.arg(filter)
  if(filter == "controls") {
    which.controls <- if(is.factor(x@ped$pheno)) x@ped$pheno == 1 else x@ped$pheno == 0
    st <- .Call('gg_geno_stats_snps', PACKAGE = "gaston", x@bed, rep(TRUE, ncol(x)), which.controls)$snps
    p <- (st$N0.f + 0.5*st$N1.f)/(st$N0.f + st$N1.f + st$N2.f)
    maf <- pmin(p, 1-p)
    w <- (maf < maf.threshold)
  } else {
    # filter = any
    w <- rep(FALSE, ncol(x))
    for(i in unique(x@ped$pheno)) {
      which.c <- (x@ped$pheno == i)
      st <- .Call('gg_geno_stats_snps', PACKAGE = "gaston", x@bed, rep(TRUE, ncol(x)), which.c)$snps
      p <- (st$N0.f + 0.5*st$N1.f)/(st$N0.f + st$N1.f + st$N2.f)
      maf <- pmin(p, 1-p)
      w <- w | (maf < maf.threshold)
    }
  }
  x <- select.snps(x, w & x@snps$maf > 0)
  if(is.factor(x@snps$genomic.region)) {
    if(!missing(min.nb.snps)) {
      nb.snps <- table(x@snps$genomic.region)
      keep <- names(nb.snps)[nb.snps >= min.nb.snps]
      x <- select.snps(x, x@snps$genomic.region %in% keep)
    }
    x@snps$genomic.region <- droplevels(x@snps$genomic.region)
  }
  x
}

Power <- function(alpha = 0.05, filter = c("controls", "any"), maf.threshold = 0.01, WSS = TRUE, C.ALPHA = TRUE, Beta.M = TRUE, ...) {
  x <- fff(...)
  x <- filter.rare.variants(x, filter, maf.threshold)
  pheno.pooled <- ifelse(x@ped$pheno==0, 0, 1)

  if(WSS){
    power.wss <- mean( WSS(x)$p.value < alpha ) 
    power.pooled.wss <- mean( WSS(x, group=pheno.pooled)$p.value < alpha )
  }
  else{  
    power.wss <- NA
    power.pooled.wss <- NA
  }

  if(C.ALPHA){ 
    power.calpha <- mean( C.ALPHA(x, target=50, B.max=1000)$p.value < alpha ) 
    power.pooled.calpha <- mean( C.ALPHA(x, group=pheno.pooled, target=50, B.max=1000)$p.value < alpha )
  }
  else{  
    power.calpha <- NA
    power.pooled.calpha <- NA
  }
  
  if(Beta.M){
    power.betam <- mean( Beta.M(x, target = 50, B.max=1000)$p.value < alpha ) 
    power.pooled.betam <- mean( Beta.M(x, group=pheno.pooled, target=50, B.max=1000)$p.value < alpha )
  }
  else{  
    power.betam <- NA
    power.pooled.betam <- NA
  }

  c("Power.WSS" = power.wss, "Power.pooled.WSS" = power.pooled.wss, "Power.C.ALPHA" = power.calpha, "Power.pooled.C.ALPHA" = power.pooled.calpha, "Power.Beta.M" = power.betam, "Power.pooled.Beta.M" = power.pooled.betam)
}
