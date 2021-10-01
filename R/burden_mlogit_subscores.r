burden.mlogit.subscores <- function(x, NullObject, genomic.region = x@snps$genomic.region, SubRegion=x@snps$SubRegion, burden.function = WSS, maf.threshold = 0.5, get.effect.size = FALSE, alpha = 0.05, cores = 10){
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
 
  alt.levels <- levels(NullObject$group)[levels(NullObject$group) != NullObject$ref.level]
  data.reg <- cbind(NullObject$data, score) ; rownames(data.reg) <- NULL
  data.reg <- dfidx(data.reg, varying = NULL, shape = "wide", choice = "ind.pheno")
  
  #Call burden on each large genomic.region
  RR <- mclapply(levels(genomic.region), function(reg) run.mlogit.subscores.withNull(reg = reg, pheno = NullObject$group, score = score[,grepl(colnames(score), pattern=reg), drop=FALSE], ref.level = NullObject$ref.level, alt.levels = alt.levels, covar.toinclude = NullObject$covar.toinclude, data = data.reg, alpha = alpha, H0.LogLik = NullObject$H0.LogLik, get.effect.size = get.effect.size), mc.cores = cores)
  if(get.effect.size){
    RR.res <- do.call(rbind, lapply(RR, function(z) z$res))
    RR.res <- as.data.frame(RR.res) ;  colnames(RR.res) <- c("p.value", "n_subscores", "is.err") ; rownames(RR.res) <- levels(genomic.region)
    RR.OR <- lapply(RR, function(z){ OR.tmp <- do.call(rbind, z$OR) ; OR.tmp <- t(OR.tmp); colnames(OR.tmp) <- c("OR", "l.lower", "l.upper") ; return(OR.tmp)}) 
    names(RR.OR) <- levels(genomic.region)
    R <- list(Asso = RR.res, effect = RR.OR)
  }else{
    RR.res <- do.call(rbind, RR) ; 
    RR.res <- as.data.frame(RR.res) ; colnames(RR.res) <- c("p.value", "n_subscores", "is.err") ; rownames(RR.res) <- levels(genomic.region)
    R <- RR.res
  }

  return(R)
}  

run.mlogit.subscores.withNull <- function (reg, pheno, score, ref.level, alt.levels, covar.toinclude, data, alpha, H0.LogLik, get.effect.size){
    region <- paste(colnames(score), collapse="+")
    if (is.null(covar.toinclude)) {
      my.formula <- Formula(as.formula(paste("ind.pheno ~ 0 |", region)))
    }else {
      my.formula <- Formula(as.formula(paste("ind.pheno ~ 0 |", region, " + ", covar.toinclude)))
    }
    
    fit <- tryCatch(mlogit(my.formula, data = data, ref.level = ref.level), error = identity, warning = identity)

    if (is(fit, "error")) {
        pval <- NA
        is.err <- 1
        OR.values <- data.frame(Estimate = rep(NA, (nlevels(pheno) - 1)*ncol(score)), sd = rep(NA,(nlevels(pheno) - 1)*ncol(score)))
    }
    else {
      if (is.null(covar.toinclude)) {
        my.model <- summary(fit)
        pval <- as.numeric(my.model$lratio$p.value)
        if (get.effect.size == TRUE) OR.values <- my.model$CoefTable
      }else{
        my.model.H1 <- summary(fit)
        pval <- pchisq(-2 * H0.LogLik + 2 * as.numeric(my.model.H1$logLik), (nlevels(pheno) - 1) * ncol(score), lower.tail = FALSE)
        if (get.effect.size == TRUE) OR.values <- my.model.H1$CoefTable
      }
      is.err <- 0
    }
    quantile.alpha <- qnorm(alpha/2, lower.tail = FALSE)
    if (get.effect.size){
      OR.values.estimate <- OR.values[grep(rownames(OR.values), pattern = reg), 1]
      OR.values.sd <- OR.values[grep(rownames(OR.values), pattern = reg), 2]
      results <- list(res=c(pval, ncol(score), is.err), OR = list(exp(OR.values.estimate), exp(OR.values.estimate - quantile.alpha * OR.values.sd), exp(OR.values.estimate + quantile.alpha * OR.values.sd)))
    }else{ results <- c(pval,  ncol(score), is.err)}
    
    #Cleaning temporary objects
    rm(score) ; rm(data) ; rm(fit) ; gc()
    return(results)
}

