burden.continuous.subscores <- function(x, NullObject, genomic.region = x@snps$genomic.region, SubRegion=x@snps$SubRegion, burden.function = WSS, maf.threshold = 0.5, get.effect.size = FALSE, alpha = 0.05, cores = 10){
  
  if (!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)
  if (missing(x)) stop("a bed.matrix 'x' is needed to compute the score")
  if(is.null(SubRegion)) warning("'SubRegion' is empty, no subscores will be included in the analysis, use rather burden()")
  
  #Temporary regions with names of big and sub regions
  x@snps$Regiontmp <- paste(genomic.region, SubRegion, sep = ".") ; x@snps$Regiontmp <- factor(x@snps$Regiontmp, levels = unique(x@snps$Regiontmp))
  
  if(!is.function(burden.function)) stop("'burden.function' should be a function. Use for example 'WSS' or 'CAST'")
  score <- burden.function(x, genomic.region = x@snps$Regiontmp)
  
  score <- as.data.frame(score)
  old.names <- colnames(score)
  names(score) <- make.names(names(score))

  # preparation data / formula
  data.reg <- cbind(NullObject$data, score) ; rownames(data.reg) <- NULL

  RR <- mclapply( levels(genomic.region), function(reg) run.burden.continuous.subscores(reg = reg, pheno = NullObject$pheno, score = score[,grepl(colnames(score), pattern=reg), drop=FALSE], covar.toinclude = NullObject$covar.toinclude, data = data.reg, get.effect.size = get.effect.size, alpha = alpha), mc.cores = cores)

  if(get.effect.size){
    RR.res <- do.call(rbind, lapply(RR, function(z) z$res))
    RR.res <- as.data.frame(RR.res) ;  colnames(RR.res) <- c("p.value", "n_subscores", "is.err") ; rownames(RR.res) <- levels(genomic.region)
    RR.beta <- lapply(RR, function(z){ beta.tmp <- do.call(rbind, z$beta) ; beta.tmp <- t(beta.tmp); colnames(beta.tmp) <- c("beta", "l.lower", "l.upper") ; return(beta.tmp)}) 
    names(RR.beta) <- levels(genomic.region)
    R <- list(Asso = RR.res, effect = RR.beta)
  }else{
    RR.res <- do.call(rbind, RR) ; 
    RR.res <- as.data.frame(RR.res) ; colnames(RR.res) <- c("p.value", "n_subscores", "is.err") ; rownames(RR.res) <- levels(genomic.region)
    R <- RR.res
  }

  return(R)

}


run.burden.continuous.subscores <- function(reg, pheno, score, region, covar.toinclude, data, get.effect.size, alpha){
  region <- paste(colnames(score), collapse="+")
  # Formula for the current region
  if(is.null(covar.toinclude)) { 
    my.formula <- as.formula(paste("ind.pheno ~ ", region))
  } else {
    my.formula <- as.formula( paste("ind.pheno ~ ", region, " + ", covar.toinclude )) 
  }

 
  # Catch errors
  fit <- tryCatch(lm(my.formula, data = data), error = identity, warning = identity)

  if(is(fit, "error")) {
    pval <- NA ; 
    is.err <- 1 ; 
    beta.values <- data.frame(Estimate = NA, sd = NA)
  } else {
    my.model <- summary(fit)
    #p-valeur totale du modele
    pval <- with(my.model, pf(fstatistic[1],fstatistic[2],fstatistic[3],lower.tail=F))
    is.err <- 0
    if(get.effect.size)  beta.values <- my.model$coefficients
  }
  
  quantile.alpha <- qnorm(alpha/2, lower.tail = FALSE)
  if (get.effect.size){
    beta.values.estimate <- beta.values[grep(rownames(beta.values), pattern = reg),2]
    beta.values.sd <- beta.values[grep(rownames(beta.values), pattern = reg), 2]
      results <- list(res=c(pval, ncol(score), is.err), beta = list(beta.values.estimate, beta.values.estimate - quantile.alpha * beta.values.sd, beta.values.estimate + quantile.alpha * beta.values.sd))
    }else{ results <- c(pval,  ncol(score), is.err)}
    
    #Cleaning temporary objects
    rm(score) ; rm(data) ; rm(fit) ; gc()
    return(results)

}
  

  
