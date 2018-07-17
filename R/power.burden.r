power.burden <- function(alpha = 2.5e-6, filter=c("whole", "controls", "any"), maf.threshold = 0.01, 
              CAST = TRUE, WSS = TRUE, pooled.analysis=TRUE, non.pooled.analysis=TRUE, analysis.by.group=TRUE,
              file.pop.maf=Kryukov, size=c(1000, 500, 500), baseline=c(0.001, 0.001), replicates=1000, select.gene=NULL, 
              same.variant=FALSE, genetic.model="additive",
              GRR.matrix, prop.del=0.5, prop.pro=0,              
              covariates=NULL, OR.alpha=NULL, reflevel="0", get.OR.value=FALSE){

##Select MAF from the file given by the user  
  if(nlevels(file.pop.maf$gene)>1){
    pop.maf <- subset(file.pop.maf, file.pop.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- file.pop.maf$maf
  }
  
#Test on GRR matrix
  if(!is.list(GRR.matrix)){
    if(is.matrix(GRR.matrix)){
      GRR.matrix <- list(GRR.matrix)
    }else{
      stop("GRR.matrix should be a list or a matrix")
    }
  }
  if(is.list(GRR.matrix)){
  #Error if general model but only one value for GRR
    if(length(GRR.matrix)==1){
      if(genetic.model=="general"){
        stop("Needs two GRR matrices in the general model")
      }else{
        GRR <- GRR.matrix[[1]]
        GRR.2 <- NULL
      }
    }else{
      if(genetic.model=="general"){
        GRR <- GRR.matrix[[1]]
        GRR.2 <- GRR.matrix[[2]]
      }else{
        warning("Only one GRR matrix needed for this model, only the first one is used")
        GRR <- GRR.matrix[[1]]
        GRR.2 <- NULL
      }
    }
  }

  
##Order arguments
  GRR.1.pars <- list(OR.del=GRR, OR.pro=1/GRR, prob.del=prop.del, prob.pro=prop.pro)
  if(!is.null(GRR.2)){
    GRR.2.pars <- list(OR.del=GRR.2, OR.pro=1/GRR.2, prob.del=prop.del, prob.pro=prop.pro)
  }else{
    GRR.2.pars <- NULL
  }
##Si meme variants: prendre fonction correspondante
  if(same.variant==FALSE){
    variant.function <- oz:::OR.matrix.fix
  }else{
    variant.function <- oz:::OR.matrix.same.fix.variant
  }
  model.pars <- list(pop.maf=pop.maf, size=size, baseline=baseline, replicates=replicates, GRR.pars=GRR.1.pars, GRR.pars.2=GRR.2.pars, OR.function=variant.function, model=genetic.model)
  
##Simulations des donnees
  x <- do.call(random.bed.matrix.GRR, model.pars)
  x <- filter.rare.variants(x, filter=filter, maf.threshold=maf.threshold)
  pheno.pooled <- ifelse(x@ped$pheno==0, 0, 1)
   	
  if(CAST){
    if(non.pooled.analysis){
      CAST.pval <- score.reg.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "CAST", alpha=OR.alpha, covariates=covariates, maf.threshold=maf.threshold, get.OR.value = get.OR.value)
      power.CAST <- mean(CAST.pval[CAST.pval$is.err == 0, "p.value"] < alpha)
      se.CAST <- sqrt((power.CAST * (1-power.CAST)) / nrow(CAST.pval[CAST.pval$is.err == 0,]))
    }else{
      power.CAST <- se.CAST <- NA
    }
    if(pooled.analysis){
      pooled.CAST.pval <- score.reg.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "CAST", alpha=OR.alpha, covariates=covariates, maf.threshold=maf.threshold, get.OR.value = get.OR.value)
      power.pooled.CAST <- mean(pooled.CAST.pval[pooled.CAST.pval$is.err==0, "p.value"] < alpha)
      se.pooled.CAST <- sqrt((power.pooled.CAST * (1-power.pooled.CAST)) / nrow(pooled.CAST.pval[pooled.CAST.pval$is.err == 0,]))
    }else{
      power.pooled.CAST <- se.pooled.CAST <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.CAST.pval <- lapply(cases, function(z) score.reg.mlogit(z, burden.score="CAST", reflevel="0"))
      power.cases.CAST <- unlist(lapply(cases.CAST.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.CAST <- sqrt(power.cases.CAST * (1-power.cases.CAST)) / nrow(power.cases.CAST)
    }else{
      power.cases.CAST <- se.cases.CAST <- rep(NA, length(baseline))
    }
  }
   
  
  if(WSS){
    if(non.pooled.analysis){
      WSS.pval <- score.reg.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "WSS", alpha=OR.alpha, covariates=covariates, maf.threshold=maf.threshold, get.OR.value = get.OR.value)
      power.WSS <- mean(WSS.pval[WSS.pval$is.err == 0, "p.value"] < alpha)
      se.WSS <- sqrt((power.WSS * (1-power.WSS)) / nrow(WSS.pval[WSS.pval$is.err == 0,]))
    }else{
      power.WSS <- se.WSS <- NA
    }
    if(pooled.analysis){
      pooled.WSS.pval <- score.reg.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, reflevel = reflevel, burden.score = "WSS", alpha=OR.alpha, covariates=covariates, maf.threshold=maf.threshold, get.OR.value = get.OR.value)
      power.pooled.WSS <- mean(pooled.WSS.pval[pooled.WSS.pval$is.err==0, "p.value"] < alpha)
      se.pooled.WSS <- sqrt((power.pooled.WSS * (1-power.pooled.WSS)) / nrow(pooled.WSS.pval[pooled.WSS.pval$is.err == 0,]))
    }else{
      power.pooled.WSS <- se.pooled.WSS <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.WSS.pval <- lapply(cases, function(z) score.reg.mlogit(z, burden.score="WSS", reflevel="0"))
      power.cases.WSS <- unlist(lapply(cases.WSS.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.WSS <- sqrt(power.cases.WSS * (1-power.cases.WSS)) / nrow(power.cases.WSS)
    }else{
      power.cases.WSS <- se.cases.WSS <- rep(NA, length(baseline))
    }
  }
  
  data.frame("power" = c(power.CAST, power.cases.CAST, power.pooled.CAST, power.WSS, power.cases.WSS, power.pooled.WSS),
  			 "se" = c(se.CAST, se.cases.CAST, se.pooled.CAST, se.WSS, se.cases.WSS, se.pooled.WSS),
  			 #Nb replicates: NA si analyse non demandee, sinon nombre de simus ayant converge
         nb.replicates = c(ifelse(CAST & non.pooled.analysis, nrow(CAST.pval[CAST.pval$is.err == 0,]), NA), 
                           if(CAST & analysis.by.group) unlist(lapply(cases.CAST.pval, function(z) nrow(z[z$is.err==0,]))) else power.cases.CAST,
                           ifelse(CAST & pooled.analysis, nrow(pooled.CAST.pval[pooled.CAST.pval$is.err == 0,]), NA),
                           ifelse(WSS & non.pooled.analysis, nrow(WSS.pval[WSS.pval$is.err == 0,]), NA),
                           if(WSS & analysis.by.group) unlist(lapply(cases.WSS.pval, function(z) nrow(z[z$is.err==0,]))) else power.cases.WSS,
                           ifelse(WSS & pooled.analysis, nrow(pooled.WSS.pval[pooled.WSS.pval$is.err == 0,]), NA)),
  			 row.names = c("CAST", paste(paste("Cases", 1:length(baseline), sep=""), "vsControls.CAST", sep=""), "pooled.CAST", "WSS", paste(paste("Cases", 1:length(baseline), sep=""), "vsControls.WSS", sep=""), "pooled.WSS")
 			)
}
		
   
