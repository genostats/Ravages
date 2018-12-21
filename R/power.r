power <- function(alpha = 2.5e-6, filter=c("whole", "controls", "any"), min.nb.snps, maf.threshold = 0.01, 
                  CAST = TRUE, WSS = TRUE, other.score, 
                  pooled.analysis=TRUE, non.pooled.analysis=TRUE, analysis.by.group=TRUE,
                  formula=NULL, data=NULL, ref.level="0",
                  model.pars=list(genes.maf=Kryukov, select.gene="R1", size=c(1000, 500, 500),
                                  baseline=c(0.001, 0.001), n.case.groups=2, replicates=1000,
                                  GRR="SKAT", GRR.multiplicative.factor=2,
                                  same.variant=FALSE, fixed.variant.prop = TRUE,
                                  genetic.model="multiplicative",
                                  prop.del=0.5, prop.pro=0)){

# dans model.pars 
#              genes.maf=Kryukov, select.gene, 
#              size=c(1000, 500, 500), baseline=c(0.001, 0.001), replicates=1000, 
#              same.variant=FALSE, fixed.variant.prop = TRUE,
#              genetic.model=c("multiplicative", "general", "recessive", "dominant"),
#              GRR.matrix.del, GRR.matrix.pro, prop.del=0.5, prop.pro=0,              
# pour compute GRR
# n.case.groups = length(baseline) (si fournis)
# GGR, GRR.value, GRR.function, GRR.mutiplicative.factor, select.gene

##Check if another score is asked
  if(!missing(other.score)){
    if(!is.function(other.score)) stop("other.score needs to be a function depending on a bed.matrix")
  }      
  
##Data simulation

  GRRmat.args <- model.pars[ intersect( names(model.pars), names(formals(GRR.matrix))) ]
  random.bm.args <- model.pars[ intersect( names(model.pars), names(formals(random.bed.matrix))) ]
  if( !("GRR.matrix.del" %in% names(random.bm.args) )) {
    GRRmat <- do.call( GRR.matrix, GRRmat.args)
    random.bm.args$GRR.matrix.del <- GRRmat
  }

  x <- do.call(random.bed.matrix, random.bm.args)

  x <- filter.rare.variants(x, filter=filter, maf.threshold=maf.threshold, min.nb.snps = min.nb.snps)
  pheno.pooled <- ifelse(x@ped$pheno==0, 0, 1)
  if (!missing(other.score)) score <- other.score(x)
   	
  if(CAST){
    if(non.pooled.analysis){
      CAST.pval <- burden.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, ref.level = ref.level, burden = "CAST", formula=formula, data=data, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.CAST <- mean(CAST.pval[CAST.pval$is.err == 0, "p.value"] < alpha)
      se.CAST <- sqrt((power.CAST * (1-power.CAST)) / nrow(CAST.pval[CAST.pval$is.err == 0,]))
    }else{
      power.CAST <- se.CAST <- NA
    }
    if(pooled.analysis){
      pooled.CAST.pval <- burden.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, ref.level = ref.level, burden = "CAST", formula=formula, data=data, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.pooled.CAST <- mean(pooled.CAST.pval[pooled.CAST.pval$is.err==0, "p.value"] < alpha)
      se.pooled.CAST <- sqrt((power.pooled.CAST * (1-power.pooled.CAST)) / nrow(pooled.CAST.pval[pooled.CAST.pval$is.err == 0,]))
    }else{
      power.pooled.CAST <- se.pooled.CAST <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(model.pars$baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.CAST.pval <- lapply(cases, function(z) burden.mlogit(z, burden="CAST", ref.level="0", formula=formula, data=data, maf.threshold=maf.threshold))
      power.cases.CAST <- unlist(lapply(cases.CAST.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.CAST <- unlist(lapply(cases.CAST.pval, function(z) sqrt(( mean(z[z$is.err==0, "p.value"]<alpha) * (1-mean(z[z$is.err==0, "p.value"]<alpha)) ) / nrow(z[z$is.err==0,]))))
    }else{
      power.cases.CAST <- se.cases.CAST <- rep(NA, length(model.pars$baseline))
    }    
  }else{
    power.CAST <- se.CAST <- power.pooled.CAST <- se.pooled.CAST <- NA
    power.cases.CAST <- se.cases.CAST <- rep(NA, length(model.pars$baseline))
  }
   
  
  if(WSS){
    if(non.pooled.analysis){
      WSS.pval <- burden.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, ref.level = ref.level, burden = "WSS", formula=formula, data=data, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.WSS <- mean(WSS.pval[WSS.pval$is.err == 0, "p.value"] < alpha)
      se.WSS <- sqrt((power.WSS * (1-power.WSS)) / nrow(WSS.pval[WSS.pval$is.err == 0,]))
    }else{
      power.WSS <- se.WSS <- NA
    }
    if(pooled.analysis){
      pooled.WSS.pval <- burden.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, ref.level = ref.level, burden = "WSS", formula=formula, data=data, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.pooled.WSS <- mean(pooled.WSS.pval[pooled.WSS.pval$is.err==0, "p.value"] < alpha)
      se.pooled.WSS <- sqrt((power.pooled.WSS * (1-power.pooled.WSS)) / nrow(pooled.WSS.pval[pooled.WSS.pval$is.err == 0,]))
    }else{
      power.pooled.WSS <- se.pooled.WSS <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(model.pars$baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.WSS.pval <- lapply(cases, function(z) burden.mlogit(z, burden="WSS", ref.level="0", formula=formula, data=data, maf.threshold=maf.threshold))
      power.cases.WSS <- unlist(lapply(cases.WSS.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.WSS <- unlist(lapply(cases.WSS.pval, function(z) sqrt(( mean(z[z$is.err==0, "p.value"]<alpha) * (1-mean(z[z$is.err==0, "p.value"]<alpha)) ) / nrow(z[z$is.err==0,]))))
    }else{
      power.cases.WSS <- se.cases.WSS <- rep(NA, length(model.pars$baseline))
    }
  }else{
    power.WSS <- se.WSS <- power.pooled.WSS <- se.pooled.WSS <- NA
    power.cases.WSS <- se.cases.WSS <- rep(NA, length(model.pars$baseline))
  }  
  
  if(!missing(other.score)){
    burden <- other.score(x) ; colnames(burden) <- levels(x@snps$genomic.region) ; rownames(burden) <- x@ped$id
    if(non.pooled.analysis){
      other.score.pval <- burden.mlogit(x, group=x@ped$pheno, genomic.region=x@snps$genomic.region, ref.level = ref.level, burden = burden, formula=formula, data=data, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.other.score <- mean(other.score.pval[other.score.pval$is.err == 0, "p.value"] < alpha)
      se.other.score <- sqrt((power.other.score * (1-power.other.score)) / nrow(other.score.pval[other.score.pval$is.err == 0,]))
    }else{
      power.other.score <- se.other.score <- NA
    }
    if(pooled.analysis){
      pooled.other.score.pval <- burden.mlogit(x, group=pheno.pooled, genomic.region=x@snps$genomic.region, ref.level = ref.level, burden = burden, formula=formula, data=data, maf.threshold=maf.threshold, get.OR.value = FALSE)
      power.pooled.other.score <- mean(pooled.other.score.pval[pooled.other.score.pval$is.err == 0, "p.value"] < alpha)
      se.pooled.other.score <- sqrt((power.pooled.other.score * (1-power.pooled.other.score)) / nrow(pooled.other.score.pval[pooled.other.score.pval$is.err == 0,]))
    }else{
      power.pooled.other.score <- se.pooled.other.score <- NA
    }
    if(analysis.by.group){
      cases <- lapply(1:length(model.pars$baseline), function(z) select.inds(x, x@ped$pheno %in% c(0,z)))
      cases.other.score.pval <- lapply(cases, function(z) burden.mlogit(z, burden=burden[z@ped$id,], ref.level="0", formula=formula, data=data, maf.threshold=maf.threshold))
      power.cases.other.score <- unlist(lapply(cases.other.score.pval, function(z) mean(z[z$is.err==0, "p.value"]<alpha)))
      se.cases.other.score <- unlist(lapply(cases.other.score.pval, function(z) sqrt(( mean(z[z$is.err==0, "p.value"]<alpha) * (1-mean(z[z$is.err==0, "p.value"]<alpha)) ) / nrow(z[z$is.err==0,]))))
    }else{
      power.cases.other.score <- se.cases.other.score <- rep(NA, length(model.pars$baseline))
    }
  }else{
    power.other.score <- se.other.score <- power.pooled.other.score <- se.pooled.other.score <- NA
    power.cases.other.score <- se.cases.other.score <- rep(NA, length(model.pars$baseline))
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
                           ifelse(!missing(other.score) & non.pooled.analysis, nrow(other.score.pval[other.score.pval$is.err == 0,]), NA),
                           if(!missing(other.score) & analysis.by.group) unlist(lapply(cases.other.score.pval, function(z) nrow(z[z$is.err==0,]))) else power.cases.other.score,
                           ifelse(!missing(other.score) & pooled.analysis, nrow(pooled.other.score.pval[pooled.other.score.pval$is.err == 0,]), NA)),
  			 row.names = c("CAST", paste(paste("Cases", 1:length(model.pars$baseline), sep=""), "vsControls.CAST", sep=""), "pooled.CAST", 
  			 			   "WSS", paste(paste("Cases", 1:length(model.pars$baseline), sep=""), "vsControls.WSS", sep=""), "pooled.WSS",
  			 			   "other.score", paste(paste("Cases", 1:length(model.pars$baseline), sep=""), "vsControls.other.score", sep=""), "pooled.other.score")
 			)
}
		
   
