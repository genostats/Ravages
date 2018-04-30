
# il faut ajouter un argument pour la fonction à appeler pour générer les OR
# (ici OR.matrix...)
random.bed.matrix <- function(pop.maf, size, baseline, replicates, OR.pars, OR.function) {
  OR.pars$n.variants <- length(pop.maf)
  nb_snps <- OR.pars$n.variants * replicates
  nb_inds <- sum(size)
  x <- new.bed.matrix(nb_inds, nb_snps);
  for(b in 1:replicates) {
    OR <- do.call( OR.function, OR.pars)
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
Power <- function(alpha = 0.05, filter = c("whole", "controls", "any"), 
                  maf.threshold = 0.01, CAST=TRUE, WSS = TRUE, C.ALPHA = TRUE, Beta.M = TRUE, 
                  Beta.M.rect = TRUE, Beta.M.freq = TRUE, SKAT=TRUE, FST = TRUE, model.pars) {
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
  
  if(FST){
    power.fst <- mean( Sum.Fst(x, target = 50, B.max = 50/alpha)$p.value < alpha ) 
    se.fst <- sqrt((power.fst * (1-power.fst))/nlevels(x@snps$genomic.region))
    power.pooled.fst <- mean( Sum.Fst(x, group=pheno.pooled, target=50, B.max = 50/alpha)$p.value < alpha )
    se.pooled.fst <- sqrt((power.pooled.fst * (1-power.pooled.fst))/nlevels(x@snps$genomic.region))
  } else {  
    power.fst <- NA ; se.fst <- NA
    power.pooled.fst <- NA ; se.pooled.fst <- NA
  }
   
  if(Beta.M){
    power.betam <- mean( Beta.M(x, target = 50, B.max = 50/alpha)$p.value < alpha ) 
    se.betam <- sqrt((power.betam * (1-power.betam))/nlevels(x@snps$genomic.region))
    power.pooled.betam <- mean( Beta.M(x, group=pheno.pooled, target=50, B.max = 50/alpha)$p.value < alpha )
    se.pooled.betam <- sqrt((power.pooled.betam * (1-power.pooled.betam))/nlevels(x@snps$genomic.region))
  }
  else{  
    power.betam <- NA ; se.betam <- NA
    power.pooled.betam <- NA ; se.pooled.betam <- NA
  }
  
  if(Beta.M.rect){
    power.betam.rect <- mean( Beta.M.rect(x, target = 50, B.max = 50/alpha)$p.value < alpha ) 
    se.betam.rect <- sqrt((power.betam.rect * (1-power.betam.rect))/nlevels(x@snps$genomic.region))
    power.pooled.betam.rect <- mean( Beta.M.rect(x, group=pheno.pooled, target=50, B.max = 50/alpha)$p.value < alpha )
    se.pooled.betam.rect <- sqrt((power.pooled.betam.rect * (1-power.pooled.betam.rect))/nlevels(x@snps$genomic.region))
  }
  else{  
    power.betam.rect <- NA ; se.betam.rect <- NA
    power.pooled.betam.rect <- NA ; se.pooled.betam.rect <- NA
  }
  
  if(Beta.M.freq){
    power.betam.freq <- mean( Beta.M.freq(x, target = 50, B.max = 50/alpha)$p.value < alpha ) 
    se.betam.freq <- sqrt((power.betam.freq * (1-power.betam.freq))/nlevels(x@snps$genomic.region))
    power.pooled.betam.freq <- mean( Beta.M.freq(x, group=pheno.pooled, target=50, B.max = 50/alpha)$p.value < alpha )
    se.pooled.betam.freq <- sqrt((power.pooled.betam.freq * (1-power.pooled.betam.freq))/nlevels(x@snps$genomic.region))
  }
  else{  
    power.betam.freq <- NA ; se.betam.freq <- NA
    power.pooled.betam.freq <- NA ; se.pooled.betam.freq <- NA
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
                        power.betam, power.pooled.betam, power.betam.rect, power.pooled.betam.rect,  power.betam.freq, power.pooled.betam.freq,
                        power.skat, power.skato, power.fst, power.pooled.fst), 
                "se"= c(se.cast, se.pooled.cast, se.wss, se.pooled.wss, se.calpha, se.pooled.calpha, 
                        se.betam, se.pooled.betam, se.betam.rect, se.pooled.betam.rect, se.betam.freq, se.pooled.betam.freq,
                        se.skat, se.skato, se.fst, se.pooled.fst), 
    "nb.replicates" = nlevels(x@snps$genomic.region),
          row.names = c("CAST", "pooled.CAST", "WSS", "pooled.WSS", "C.alpha", "pooled.C.alpha", 
                        "Beta.M", "pooled.Beta.M", "Beta.M.rect", "pooled.Beta.M.rect", "Beta.M.freq", "pooled.Beta.M.freq",
                        "SKAT", "SKAT-O", "Sum.Fst", "pooled.Sum.Fst"))  

}


power.burden <- function(alpha = 0.05, filter=c("whole", "controls", "any"), maf.threshold = 0.01, CAST = TRUE, WSS = TRUE, burden = FALSE, regression = TRUE, model.pars){
  if(is.character(model.pars[[1]])){
  	model.pars[[1]] <- GnomADgenes$maf[GnomADgenes$gene==model.pars[[1]]]
  }
  x <- do.call(random.bed.matrix, model.pars)
  x <- filter.rare.variants(x, filter, maf.threshold)
  pheno.pooled <- ifelse(x@ped$pheno==0, 0, 1)
	
  if(CAST){
    if(burden){
      power.CAST.burden <- mean(CAST(x, maf.threshold = maf.threshold)$p.value < alpha)
      se.CAST.burden <- sqrt((power.CAST.burden * (1-power.CAST.burden)) / nlevels(x@snps$genomic.region))
      power.pooled.CAST.burden <- mean(CAST(x, maf.threshold = maf.threshold, group = pheno.pooled)$p.value < alpha)
      se.pooled.CAST.burden <- sqrt((power.pooled.CAST.burden * (1-power.pooled.CAST.burden)) / nlevels(x@snps$genomic.region))
    }else{
      power.CAST.burden <- se.CAST.burden <- power.pooled.CAST.burden <- se.pooled.CAST.burden <- NA
    }
    
    if(regression){
      CAST.regression <- score.reg.mlogit(x, reflevel = "1", burden.score = "CAST", get.OR.value = FALSE)
      power.CAST.regression <- mean(CAST.regression[CAST.regression$is.err == 0, "p.value"] < alpha)
      se.CAST.regression <- sqrt(power.CAST.regression * (1-power.CAST.regression)) / nrow(CAST.regression[CAST.regression$is.err == 0]))
      pooled.CAST.regression <- score.reg.mlogit(x, group = pheno.pooled, reflevel = "1", burden.score = "CAST", get.OR.value = FALSE)
      power.pooled.CAST.regresion <- mean(pooled.CAST.regression[pooled.CAST.regression$is.err==0, "p.value"] < alpha)
      se.pooled.CAST.regression <- sqrt((power.pooled.CAST.regression * (1-power.pooled.CAST.regression)) / nrow(pooled.CAST.regression[pooled.CAST.regression$is.err == 0]))
    }else{
      power.CAST.regression <- se.CAST.regression <- power.pooled.CAST.regression <- se.pooled.CAST.regression <- NA
    }
  }else{
    power.CAST.burden <- se.CAST.burden <- power.pooled.CAST.burden <- se.pooled.CAST.burden <- NA
    power.CAST.regression <- se.CAST.regression <- power.pooled.CAST.regression <- se.pooled.CAST.regression <- NA
  }
    
  
  if(WSS){
    if(burden){
      power.WSS.burden <- mean(WSS(x, maf.threshold = maf.threshold)$p.value < alpha)
      se.WSS.burden <- sqrt((power.WSS.burden * (1-power.WSS.burden)) / nlevels(x@snps$genomic.region))
      power.pooled.WSS.burden <- mean(WSS(x, maf.threshold = maf.threshold, group = pheno.pooled)$p.value < alpha)
      se.pooled.WSS.burden <- sqrt((power.pooled.WSS.burden * (1-power.pooled.WSS.burden)) / nlevels(x@snps$genomic.region))
    }else{
      power.WSS.burden <- se.WSS.burden <- power.pooled.WSS.burden <- se.pooled.WSS.burden <- NA
    }
    
    if(regression){
      WSS.regression <- score.reg.mlogit(x, reflevel = "1", burden.score = "WSS", get.OR.value = FALSE)
      power.WSS.regression <- mean(WSS.regression[WSS.regression$is.err == 0, "p.value"] < alpha)
      se.WSS.regression <- sqrt(power.WSS.regression * (1-power.WSS.regression)) / nrow(WSS.regression[WSS.regression$is.err == 0]))
      pooled.WSS.regression <- score.reg.mlogit(x, group = pheno.pooled, reflevel = "1", burden.score = "WSS", get.OR.value = FALSE)
      power.pooled.WSS.regresion <- mean(pooled.WSS.regression[pooled.WSS.regression$is.err==0, "p.value"] < alpha)
      se.pooled.WSS.regression <- sqrt((power.pooled.WSS.regression * (1-power.pooled.WSS.regression)) / nrow(pooled.WSS.regression[pooled.WSS.regression$is.err == 0]))
    }else{
      power.WSS.regression <- se.WSS.regression <- power.pooled.WSS.regression <- se.pooled.WSS.regression <- NA
    }
  }else{
    power.WSS.burden <- se.WSS.burden <- power.pooled.WSS.burden <- se.pooled.WSS.burden <- NA
    power.WSS.regression <- se.WSS.regression <- power.pooled.WSS.regression <- se.pooled.WSS.regression <- NA
  }
  
  data.frame("power" = c(power.CAST.burden, power.pooled.CAST.burden, power.CAST.regression, power.pooled.CAST.regression, 
  						 power.WSS.burden, power.pooled.WSS.burden, power.WSS.regression, power.pooled.WSS.regression),
  			 "se" = c(se.CAST.burden, se.pooled.CAST.burden, se.CAST.regression, se.pooled.CAST.regression, 
  			 		  se.WSS.burden, se.pooled.WSS.burden, se.WSS.regression, se.pooled.WSS.regression),
  			 nb.replicates = c(rep(nlevels(x@snps$genomic.region), 2), nrow(CAST.regression[CAST.regression$is.err == 0]), nrow(pooled.CAST.regression[pooled.CAST.regression$is.err == 0]),
  			 				   rep(nlevels(x@snps$genomic.region), 2), nrow(WSS.regression[WSS.regression$is.err == 0]), nrow(pooled.WSS.regression[pooled.WSS.regression$is.err == 0]) ),
  			 row.names = c("CAST.burden", "pooled.CAST.burden", "CAST.regression", "pooled.CAST.regression",
  			 			   "WSS.burden", "pooled.WSS.burden", "WSS.regression", "pooled.WSS.regression")
 			)
}
		
