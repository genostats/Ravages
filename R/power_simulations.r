
# il faut ajouter un argument pour la fonction à appeler pour générer les OR
# (ici OR.matrix...)
random.bed.matrix <- function(pop.maf, size, baseline, replicates, OR.pars, scenario=c(1,2)) {
  OR.pars$n.variants <- length(pop.maf)
  nb_snps <- OR.pars$n.variants * replicates
  nb_inds <- sum(size)
  x <- new.bed.matrix(nb_inds, nb_snps);
  for(b in 1:replicates) {
    if(scenario == 1){
      OR <- do.call( same.OR.matrix, OR.pars)
    } else {
      OR <- do.call( OR.matrix, OR.pars)
    }
    MAFS <- group.mafs(pop.maf, OR, baseline)
    .Call("oz_random_filling_bed_matrix", PACKAGE = "oz", x@bed, MAFS, size, (b-1)*OR.pars$n.variants)
  }
  x@ped$pheno <- rep.int( 1:length(size) - 1, size)
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = OR.pars$n.variants) )
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x <- set.stats(x, verbose = FALSE)
  x
}


filter.rare.variants <- function(x, filter = c("whole", "controls", "any"), maf.threshold = 0.01, min.nb.snps) {
  filter <- match.arg(filter)
  if(filter == "controls") {
    which.controls <- if(is.factor(x@ped$pheno)) x@ped$pheno == 1 else x@ped$pheno == 0
    st <- .Call('gg_geno_stats_snps', PACKAGE = "gaston", x@bed, rep(TRUE, ncol(x)), which.controls)$snps
    p <- (st$N0.f + 0.5*st$N1.f)/(st$N0.f + st$N1.f + st$N2.f)
    maf <- pmin(p, 1-p)
    w <- (maf < maf.threshold)
  }
  if(filter == "any"){
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
  if(filter == "whole"){
      ##Filter = whole
      w <- (x@snps$maf < maf.threshold)
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

# model.pars doit contenir les arguments de random.bed.matrix.with.model
Power <- function(alpha = 0.05, filter = c("whole", "controls", "any"), maf.threshold = 0.01, CAST=TRUE, WSS = TRUE, C.ALPHA = TRUE, Beta.M = TRUE, SKAT=TRUE, model.pars) {
  x <- do.call(random.bed.matrix, model.pars)
  x <- filter.rare.variants(x, filter, maf.threshold)
  pheno.pooled <- ifelse(x@ped$pheno==0, 0, 1)
  
  if(CAST){
    power.cast <- mean(CAST(x, maf.threshold = maf.threshold)$p.value < alpha)
    se.cast <- sqrt((power.cast * (1-power.cast)) / nlevels(x@snps$genomic.region))
    power.pooled.cast <- mean(CAST(x, group = pheno.pooled, maf.threshold = maf.threshold)$p.value < alpha)
    se.pooled.cast <- sqrt((power.pooled.cast * (1-power.pooled.cast)) / nlevels(x@snps$genomic.region))
  }
  else {
    power.cast <- NA ; se.cast <- NA
    power.pooled.cast <- NA ; se.pooled.cast <- NA
  }
  
  if(WSS){
    power.wss <- mean( WSS(x)$p.value < alpha ) 
    se.wss <- sqrt((power.wss*(1-power.wss))/nlevels(x@snps$genomic.region)) 
    power.pooled.wss <- mean( WSS(x, group=pheno.pooled)$p.value < alpha )
    se.pooled.wss <- sqrt((power.pooled.wss*(1-power.pooled.wss))/nlevels(x@snps$genomic.region))
  }
  else{  
    power.wss <- NA ; se.wss <- NA
    power.pooled.wss <- NA ; se.pooled.wss <- NA
  }

  if(C.ALPHA){ 
    power.calpha <- mean( C.ALPHA(x, target=50, B.max = 50/alpha)$p.value < alpha ) 
    se.calpha <- sqrt((power.calpha*(1-power.calpha))/nlevels(x@snps$genomic.region))
    power.pooled.calpha <- mean( C.ALPHA(x, group=pheno.pooled, target=50, B.max = 50/alpha)$p.value < alpha )
    se.pooled.calpha <- sqrt((power.pooled.calpha*(1-power.pooled.calpha))/nlevels(x@snps$genomic.region))
  }
  else{  
    power.calpha <- NA ; se.calpha <- NA
    power.pooled.calpha <- NA ; se.pooled.calpha <- NA
  }
  
  if(Beta.M){
    power.betam <- mean( Beta.M(x, target = 50, B.max = 50/alpha)$p.value < alpha ) 
    se.betam <- sqrt((power.betam * (1-power.betam))/nlevels(x@snps$genomic.region))
    power.pooled.betam <- mean( Beta.M(x, group=pheno.pooled, target=50, B.max = 50/alpha)$p.value < alpha )
    se.pooled.betam <- sqrt((power.pooled.betam * (1-power.pooled.betam))/nlevels(x@snps$genomic.region))

    power.betam.rect <- mean( Beta.M.rect(x, target = 50, B.max = 50/alpha)$p.value < alpha ) 
    se.betam.rect <- sqrt((power.betam.rect * (1-power.betam.rect))/nlevels(x@snps$genomic.region))
    power.pooled.betam.rect <- mean( Beta.M.rect(x, group=pheno.pooled, target=50, B.max = 50/alpha)$p.value < alpha )
    se.pooled.betam.rect <- sqrt((power.pooled.betam.rect * (1-power.pooled.betam.rect))/nlevels(x@snps$genomic.region))
  }
  else{  
    power.betam <- NA ; se.betam <- NA
    power.pooled.betam <- NA ; se.pooled.betam <- NA
    power.betam.rect <- NA ; se.betam.rect <- NA
    power.pooled.betam.rect <- NA ; se.pooled.betam.rect <- NA
  }
  
  if (SKAT){
    obj.null <- SKAT_Null_Model(pheno.pooled ~ 1, out_type="D")
    SKAT.p <- sapply(levels(x@snps$genomic.region), function(y) SKAT(gaston:::as.matrix(select.snps(x, x@snps$genomic.region==y)), obj.null , r.corr=0)$p.value)
    SKATO.p <- sapply(levels(x@snps$genomic.region), function(y) SKAT(gaston:::as.matrix(select.snps(x, x@snps$genomic.region==y)), obj.null , method="SKATO")$p.value)
    power.skat <- mean(SKAT.p<alpha)
    power.skato <- mean(SKATO.p<alpha)
    se.skat <- sqrt((power.skat * (1-power.skat))/nlevels(x@snps$genomic.region))
    se.skato <- sqrt((power.skato * (1-power.skato))/nlevels(x@snps$genomic.region))
  }
  else{
    power.skat <- NA ; se.skat <- NA
    power.skato <- NA ; se.skato <- NA
  }
  
  data.frame(
    "power" = c(power.cast, power.pooled.cast, power.wss, power.pooled.wss, power.calpha, power.pooled.calpha, 
                power.betam, power.pooled.betam, power.betam.rect, power.pooled.betam.rect, 
                power.skat, power.skato, nlevels(x@snps$genomic.region)), 
        "se"= c(se.cast, se.pooled.cast, se.wss, se.pooled.wss, se.calpha, se.pooled.calpha, 
                se.betam, se.pooled.betam, se.betam.rect, se.pooled.betam.rect, se.skat, se.skato, NA), 
  row.names= c("CAST", "pooled.CAST", "WSS", "pooled.WSS", "C.alpha", "pooled.C.alpha", "Beta.M", "pooled.Beta.M", 
               "Beta.M.rect", "pooled.Beta.M.rect", "SKAT", "SKAT-O", "nb.replicates"))  

}
