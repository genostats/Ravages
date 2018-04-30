######Utilisation de mlogit
score.reg.mlogit <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, burden.score = c("CAST", "WSS", "Other"), other.score = NULL, reflevel, covariates=NULL, alpha=0.05){
  if(burden.score == "CAST"){
    score <- oz:::CAST.0(x, genomic.region)
  }
  if(burden.score == "WSS"){
    score <- oz:::WSS.0(x, genomic.region)
  }
  if(burden.score == "Other"){
    if(is.null(other.score)){
      stop("ERROR, need to specify a score if the score is different from CAST or WSS")
      }
    if(!is.matrix(other.score) | ncol(other.score) != nlevels(x@snps$genomic.region) | nrow(other.score) != nrow(x@ped) ){
      stop("Score is not a matrix or has wrong dimensions")
      }
    score <- other.score
  }
  
  if(!is.null(covariates)){
    if(nrow(covariates)!=length(pheno)){
      stop("Covariates has wrong dimensions")
    }
  }
  
  if(!is.factor(group)){
    group <- as.factor(group)
  }
  
  alt.levels <- levels(group)[!(levels(group) %in% reflevel)]
  
  pval <- as.data.frame(t(sapply(unique(genomic.region), function(z) get.model.parameters.mlogit(pheno = group, score = score, region=z, reflevel=reflevel, alt.levels=alt.levels, covariates=covariates, alpha=alpha))))
  colnames(pval) <- c("p.value", "is.err", paste("OR", alt.levels, sep="."), paste("l.lower", alt.levels, sep="."), paste("l.upper", alt.levels, sep="."))
  rownames(pval) <- unique(genomic.region)
  return(pval)
}


get.model.parameters.mlogit <- function(pheno = group, score = score, region, reflevel, alt.levels, covariates, alpha){
  assign("last.warning", NULL, envir = baseenv())
  
  #Model
  score.pheno <- cbind(score[,region], pheno, covariates)
  score.pheno <- as.data.frame(score.pheno)
  ifelse(is.null(covariates), colnames(score.pheno) <- c("region", "ind.pheno"), colnames(score.pheno) <- c("region", "ind.pheno", paste("covar", 1:ncol(covariates), sep="")))
  rownames(score.pheno) <- 1:nrow(score.pheno)
  score.mlogit <- mlogit.data(score.pheno, varying=NULL, shape="wide", choice="ind.pheno", alt.levels=levels(pheno))
  
  #Formula for the model
  if(is.null(covariates)){ 
    my.formula <- mFormula(ind.pheno ~ 0 | region)
  }
  else{
    if(ncol(covariates)==1){ my.formula <- mFormula(ind.pheno ~ 0 | region + covar1) ; my.formula.H0 <- mFormula(ind.pheno ~ 0 | covar1) }
    if(ncol(covariates)==2){ my.formula <- mFormula(ind.pheno ~ 0 | region + covar1 + covar2) ; my.formula.H0 <- mFormula(ind.pheno ~ 0 | covar1 + covar2) }
    if(ncol(covariates)==3){ my.formula <- mFormula(ind.pheno ~ 0 | region + covar1 + covar2 + covar3) ; my.formula.H0 <- mFormula(ind.pheno ~ 0 | covar1 + covar2 + covar3) }
    if(ncol(covariates)==4){ my.formula <- mFormula(ind.pheno ~ 0 | region + covar1 + covar2 + covar3 + covar4) ; my.formula.H0 <- mFormula(ind.pheno ~ 0 | covar1 + covar2 + covar3 + covar4) }
    if(ncol(covariates)==5){ my.formula <- mFormula(ind.pheno ~ 0 | region + covar1 + covar2 + covar3 + covar4 + covar5) ; mFormula.H0(ind.pheno ~ 0 | covar1 + covar2 + covar3 + covar4 + covar5) }
  }
  
  #Test if Error
  if(is(tryCatch(mlogit(my.formula, data=score.mlogit, reflevel = reflevel), error=function(w) w, warning=function(y) y), "error")){
    pval <- NA ; is.err <- 1
  }
  else{
    if(is.null(covariates)){
      my.model <- summary(mlogit(my.formula, data=score.mlogit, reflevel = reflevel))
      pval <- as.numeric(my.model$lratio$p.value)
      OR.values <- my.model$CoefTable
    }
    else{
      my.model.H0 <- summary(mlogit(my.formula.H0, data=score.mlogit, reflevel = reflevel))
      my.model.H1 <- summary(mlogit(my.formula, data=score.mlogit, reflevel = reflevel))
      pval <- pchisq(-2*as.numeric(my.model.H0$logLik) + 2*as.numeric(my.model.H1$logLik), nlevels(pheno)-1, lower.tail=FALSE)
      OR.values <- my.model.H1$CoefTable
    }
    is.err <- 0
  }

  quantile.alpha <- qnorm(alpha/2, lower.tail=FALSE)

  return(c(pval, is.err, 
        as.numeric(exp(OR.values[paste(alt.levels, "region", sep=":"),1])), 
        as.numeric(exp(OR.values[paste(alt.levels, "region", sep=":"),1]-quantile.alpha*OR.values[paste(alt.levels, "region", sep=":"),2])), 
        as.numeric(exp(OR.values[paste(alt.levels, "region", sep=":"),1]+quantile.alpha*OR.values[paste(alt.levels, "region", sep=":"),2]))
        ))
}
  
