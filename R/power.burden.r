power.burden <- function(alpha = 2.5e-6, filter=c("whole", "controls", "any"), min.nb.snps = NULL, maf.threshold = 0.01, 
              CAST = c(TRUE, FALSE), WSS = c(TRUE, FALSE), other.score=NULL, 
              pooled.analysis=c(TRUE, FALSE), non.pooled.analysis=c(TRUE, FALSE), analysis.by.group=c(TRUE, FALSE),
              file.pop.maf=Ravages::Kryukov, size=c(1000, 500, 500), baseline=c(0.001, 0.001), replicates=1000, select.gene=NULL, 
              same.variant=c(FALSE, TRUE), fixed.variant.prop = c(TRUE, FALSE),
              genetic.model=c("multiplicative", "general", "recessive", "dominant"),
              GRR.matrix, GRR.matrix.pro=NULL, prop.del=0.5, prop.pro=0,              
              covariates=NULL, reflevel="0"){

##Check if another score is asked
  if(!is.null(other.score)){
    if(!is.function(other.score)) stop("other.score needs to be a function depending on a bed.matrix")
  }      

##Arguments for data simulation
  model.pars <- list(file.pop.maf=file.pop.maf, size=size, baseline=baseline, replicates=replicates, GRR.matrix=GRR.matrix, GRR.matrix.pro=GRR.matrix.pro, same.variant=same.variant, fixed.variant.prop = fixed.variant.prop, genetic.model=genetic.model, select.gene=select.gene, prop.del = prop.del, prop.pro=prop.pro)
  
##Simulations des donnees
  x <- do.call(random.bed.matrix.GRR, model.pars)
  x <- filter.rare.variants(x, filter=filter, maf.threshold=maf.threshold, min.nb.snps = min.nb.snps)
  pheno.pooled <- ifelse(x@ped$pheno==0, 0, 1)
  if (!is.null(other.score)) score <- other.score(x)
   	
  if(CAST){
    if(non.pooled.analysis){
      CAST.pval <- score.reg.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "CAST", covariates=covariates, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.CAST <- mean(CAST.pval[CAST.pval$is.err == 0, "p.value"] < alpha)
      se.CAST <- sqrt((power.CAST * (1-power.CAST)) / nrow(CAST.pval[CAST.pval$is.err == 0,]))
    }else{
      power.CAST <- se.CAST <- NA
    }
    if(pooled.analysis){
      pooled.CAST.pval <- score.reg.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "CAST", covariates=covariates, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.pooled.CAST <- mean(pooled.CAST.pval[pooled.CAST.pval$is.err==0, "p.value"] < alpha)
      se.pooled.CAST <- sqrt((power.pooled.CAST * (1-power.pooled.CAST)) / nrow(pooled.CAST.pval[pooled.CAST.pval$is.err == 0,]))
    }else{
      power.pooled.CAST <- se.pooled.CAST <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.CAST.pval <- lapply(cases, function(z) score.reg.mlogit(z, burden.score="CAST", reflevel="0", maf.threshold=maf.threshold))
      power.cases.CAST <- unlist(lapply(cases.CAST.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.CAST <- unlist(lapply(cases.CAST.pval, function(z) sqrt(( mean(z[z$is.err==0, "p.value"]<alpha) * (1-mean(z[z$is.err==0, "p.value"]<alpha)) ) / nrow(z[z$is.err==0,]))))
    }else{
      power.cases.CAST <- se.cases.CAST <- rep(NA, length(baseline))
    }    
  }else{
    power.CAST <- se.CAST <- power.pooled.CAST <- se.pooled.CAST <- NA
    power.cases.CAST <- se.cases.CAST <- rep(NA, length(baseline))
  }
   
  
  if(WSS){
    if(non.pooled.analysis){
      WSS.pval <- score.reg.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "WSS", covariates=covariates, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.WSS <- mean(WSS.pval[WSS.pval$is.err == 0, "p.value"] < alpha)
      se.WSS <- sqrt((power.WSS * (1-power.WSS)) / nrow(WSS.pval[WSS.pval$is.err == 0,]))
    }else{
      power.WSS <- se.WSS <- NA
    }
    if(pooled.analysis){
      pooled.WSS.pval <- score.reg.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "WSS", covariates=covariates, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.pooled.WSS <- mean(pooled.WSS.pval[pooled.WSS.pval$is.err==0, "p.value"] < alpha)
      se.pooled.WSS <- sqrt((power.pooled.WSS * (1-power.pooled.WSS)) / nrow(pooled.WSS.pval[pooled.WSS.pval$is.err == 0,]))
    }else{
      power.pooled.WSS <- se.pooled.WSS <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.WSS.pval <- lapply(cases, function(z) score.reg.mlogit(z, burden.score="WSS", reflevel="0", maf.threshold=maf.threshold))
      power.cases.WSS <- unlist(lapply(cases.WSS.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.WSS <- unlist(lapply(cases.WSS.pval, function(z) sqrt(( mean(z[z$is.err==0, "p.value"]<alpha) * (1-mean(z[z$is.err==0, "p.value"]<alpha)) ) / nrow(z[z$is.err==0,]))))
    }else{
      power.cases.WSS <- se.cases.WSS <- rep(NA, length(baseline))
    }
  }else{
    power.WSS <- se.WSS <- power.pooled.WSS <- se.pooled.WSS <- NA
    power.cases.WSS <- se.cases.WSS <- rep(NA, length(baseline))
  }  
  
  if(!is.null(other.score)){
    if(non.pooled.analysis){
      other.score.pval <- score.reg.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "Other", other.score=score, covariates=covariates, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.other.score <- mean(other.score.pval[other.score.pval$is.err == 0, "p.value"] < alpha)
      se.other.score <- sqrt((power.other.score * (1-power.other.score)) / nrow(other.score.pval[other.score.pval$is.err == 0,]))
    }else{
      power.other.score <- se.other.score <- NA
    }
    if(pooled.analysis){
      pooled.other.score.pval <- score.reg.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "Other", other.score=score, covariates=covariates, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.pooled.other.score <- mean(pooled.other.score.pval[pooled.other.score.pval$is.err == 0, "p.value"] < alpha)
      se.pooled.other.score <- sqrt((power.pooled.other.score * (1-power.pooled.other.score)) / nrow(pooled.other.score.pval[pooled.other.score.pval$is.err == 0,]))
    }else{
      power.pooled.other.score <- se.pooled.other.score <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.other.score.pval <- lapply(cases, function(z) score.reg.mlogit(z, burden.score= "Other", other.score=other.score(z), reflevel="0", maf.threshold=maf.threshold))
      power.cases.other.score <- unlist(lapply(cases.other.score.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.other.score <- unlist(lapply(cases.other.score.pval, function(z) sqrt(( mean(z[z$is.err==0, "p.value"]<alpha) * (1-mean(z[z$is.err==0, "p.value"]<alpha)) ) / nrow(z[z$is.err==0,]))))
    }else{
      power.cases.other.score <- se.cases.other.score <- rep(NA, length(baseline))
    }
  }else{
    power.other.score <- se.other.score <- power.pooled.other.score <- se.pooled.other.score <- NA
    power.cases.other.score <- se.cases.other.score <- rep(NA, length(baseline))
  }
      
  
  data.frame("power" = c(power.CAST, power.cases.CAST, power.pooled.CAST, power.WSS, power.cases.WSS, power.pooled.WSS, power.other.score, power.cases.other.score, power.pooled.other.score),
  			 "se" = c(se.CAST, se.cases.CAST, se.pooled.CAST, se.WSS, se.cases.WSS, se.pooled.WSS, se.other.score, se.cases.other.score, se.pooled.other.score),
  			 #Nb replicates: NA si analyse non demandee, sinon nombre de simus ayant converge
         nb.replicates = c(ifelse(CAST & non.pooled.analysis, nrow(CAST.pval[CAST.pval$is.err == 0,]), NA), 
                           if(CAST & analysis.by.group) unlist(lapply(cases.CAST.pval, function(z) nrow(z[z$is.err==0,]))) else power.cases.CAST,
                           ifelse(CAST & pooled.analysis, nrow(pooled.CAST.pval[pooled.CAST.pval$is.err == 0,]), NA),
                           ifelse(WSS & non.pooled.analysis, nrow(WSS.pval[WSS.pval$is.err == 0,]), NA),
                           if(WSS & analysis.by.group) unlist(lapply(cases.WSS.pval, function(z) nrow(z[z$is.err==0,]))) else power.cases.WSS,
                           ifelse(WSS & pooled.analysis, nrow(pooled.WSS.pval[pooled.WSS.pval$is.err == 0,]), NA),
                           ifelse(!is.null(other.score) & non.pooled.analysis, nrow(other.score.pval[other.score.pval$is.err == 0,]), NA),
                           if(!is.null(other.score) & analysis.by.group) unlist(lapply(cases.other.score.pval, function(z) nrow(z[z$is.err==0,]))) else power.cases.other.score,
                           ifelse(!is.null(other.score) & pooled.analysis, nrow(pooled.other.score.pval[pooled.other.score.pval$is.err == 0,]), NA)),
  			 row.names = c("CAST", paste(paste("Cases", 1:length(baseline), sep=""), "vsControls.CAST", sep=""), "pooled.CAST", 
  			 			   "WSS", paste(paste("Cases", 1:length(baseline), sep=""), "vsControls.WSS", sep=""), "pooled.WSS",
  			 			   "other.score", paste(paste("Cases", 1:length(baseline), sep=""), "vsControls.other.score", sep=""), "pooled.other.score")
 			)
}
		
   
