
# un des scenario d'OR ...
OR.matrix <- function(n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro) {
  if(length(OR.del) != length(OR.pro)) 
    stop("Dimensions mismatch")
  OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
  # neutral, deleterious or protective
  v <- sample(1:3, n.variants, TRUE, c(1-prob.del-prob.pro, prob.del, prob.pro))
  t(apply(OR, 1, function(or) or[v]))
}
#example
#0R.matrix(20 , c(2,4), c(0.5,0.25), 0.2, 0.1)

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

filter.rare.variants <- function(x, filter = c("controls", "any"), maf.threshold = 0.01) {
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
  x <- select.snps(x, w)
  if(is.factor(x@snps$genomic.region)) 
    x@snps$genomic.region <- droplevels(x@snps$genomic.region)
  x
}

Power <- function(alpha = 0.05, filter = c("controls", "any"), maf.threshold = 0.01, WSS = TRUE, C.ALPHA = TRUE, Beta.M = TRUE, ...) {
  x <- fff(...)
  x <- filter.rare.variants(x, filter, maf.threshold)
  
  if(WSS) 
    power.wss <- mean( WSS(x)$p.value < alpha ) 
  else  
    power.wss <- NA

  if(C.ALPHA) 
    power.calpha <- mean( C.ALPHA(x, target = 10, B.max = 100)$p.value < alpha ) 
  else  
    power.calpha <- NA

  if(Beta.M) 
    power.betam <- mean( Beta.M(x, target = 10, B.max = 100)$p.value < alpha ) 
  else  
    power.betam <- NA

  c("Power.WSS" = power.wss, "Power.C.ALPHA" = power.calpha, "Power.Beta.M" = power.betam)
}
