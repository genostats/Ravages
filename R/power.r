
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
    SKAT.p <- sapply(levels(x@snps$genomic.region), function(y) SKAT(gaston::as.matrix(select.snps(x, x@snps$genomic.region==y)), obj.null , r.corr=0)$p.value)
    SKATO.p <- sapply(levels(x@snps$genomic.region), function(y) SKAT(gaston::as.matrix(select.snps(x, x@snps$genomic.region==y)), obj.null , method="SKATO")$p.value)
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

