score.reg <- function(x, pheno = x@ped$pheno, genomic.region = x@snps$genomic.region, class.ref, burden.score = c("CAST", "WSS", "Other"), other.score = NULL, covariates = NULL, print.err = TRUE){
  if(!is.factor(pheno)){
    warnings("pheno is converted to a factor")
    pheno <- as.factor(pheno)
  }
  
  if(!is.null(covariates)){
    if(nrow(covariates)!=length(pheno)){
      stop("Covariates as wrong dimensions")
    }
  }
  
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
  
  #LRT and df for each genomic region, one genomic region by line
  model.parameters <- t(apply(score, 2, function(z) get.model.parameters(pheno=pheno, score=z, class.ref = class.ref, covariates = covariates, print.err = print.err)))
  colnames(model.parameters)=c("LRT", "nb.df", "is.err")
  
  #P-value
  pval <- apply(model.parameters, 1, function(z) pchisq(z["LRT"], z["nb.df"], lower.tail=FALSE))
  
  #Results
  R <- data.frame("LRT" = model.parameters[,"LRT"], "df" = model.parameters[,"nb.df"], "p.value" = pval, is.err = model.parameters[,"is.err"])
  R
}
  
  
  
get.model.parameters <- function(pheno = pheno, score, class.ref = calss.ref, covariates = NULL, print.err){
  assign("last.warning", NULL, envir = baseenv())
  #Null model H0
  if(is.null(covariates)){
    null.model <- vglm(pheno ~ 1, family=multinomial(refLevel=class.ref))
  }
  else{
    null.model <- vglm(pheno ~ ., data=covariates, family=multinomial(refLevel=class.ref))
  }
  
  #Model to test
  if(is.null(covariates)){
    score.model <- vglm(pheno ~ score, family=multinomial(refLevel=class.ref))
  }
  else{
    score.model <- vglm(pheno ~ score + ., data=covariates, family=multinomial(refLevel=class.ref))
  }
  
  #Error
  #sans covariables
  if(print.err==TRUE){
    if(is.null(covariates)){
      is.err <- ifelse(is(tryCatch(vglm(pheno ~ score, family=multinomial(refLevel=class.ref)), warning=function(w) w), "warning"), 1, 0)
    }
    else{
      is.err <- ifelse(is(tryCatch(vglm(pheno ~ score + ., data=covariates, family=multinomial(refLevel=class.ref)), warning=function(w) w), "warning"), 1, 0)
    }
  }
  else{
    is.err <- NA
  }
  
  #LRT
  LRT <- null.model@criterion$deviance - score.model@criterion$deviance
  
  #Number of df
  nb.df <- null.model@df.residual - score.model@df.residual
  
  return(c(LRT, nb.df, is.err))
}


  