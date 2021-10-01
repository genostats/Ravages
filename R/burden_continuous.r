burden.continuous <- function(x, NullObject, genomic.region = x@snps$genomic.region, burden, maf.threshold = 0.5, get.effect.size = F, alpha = 0.05, cores = 10){
  
  if(is.numeric(burden)) {
    if(!is.matrix(burden)){
      stop("Score is not a matrix")
    } else {
      if(is.null(colnames(burden))){ 
        colnames(burden) <- make.names(1:ncol(burden))
      }
    }
    #Check between number of individuals
    if(nrow(burden) != length(NullObject$pheno)) stop("Different number of individuals in 'burden' and 'NullObject'")
    score <- burden
  } else { 
    if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
    genomic.region <- droplevels(genomic.region)
    if(missing(x)) stop("a bed.matrix 'x' is needed to compute the score")

    #Check between number of individuals
    if(nrow(x) != length(NullObject$pheno)) stop("Different number of individuals in 'x' and 'NullObject'")

    if(burden == "CAST"){
      score <- CAST(x, genomic.region, maf.threshold)
    } else if(burden == "WSS"){
      score <- WSS(x, genomic.region)
    } else {
      stop("'burden' should be \"CAST\", \"WSS\", or a matrix of pre-computed burdens");
    }
  }
  score <- as.data.frame(score)
  # to ensure syntactically correct formulas
  old.names <- colnames(score)
  names(score) <- make.names(names(score))

  # preparation data / formula
  data.reg <- cbind(NullObject$data, score) ; rownames(data.reg) <- NULL

  R <- do.call(rbind, mclapply( names(score), function(reg) run.burden.continuous(pheno = NullObject$pheno, score = score, region = reg, covar.toinclude = NullObject$covar.toinclude, data = data.reg, get.effect.size = get.effect.size, alpha = alpha), mc.cores = cores))

  R <- as.data.frame( R );

  if (get.effect.size) colnames(R) <- c("p.value", "is.err", "beta", "l.lower", "l.upper")
  else colnames(R) <- c("p.value", "is.err")

  rownames(R) <- old.names

  return(R)
}


run.burden.continuous <- function(pheno, score, region, covar.toinclude, data, get.effect.size, alpha){
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
    pval <- my.model$coefficients[region, 4]
    is.err <- 0
    if(get.effect.size)  beta.values <- my.model$coefficients
  }
  
  quantile.alpha <- qnorm(alpha/2, lower.tail = FALSE)
  if (get.effect.size){
    beta.values.estimate <- beta.values[region,2]
    beta.values.sd <- beta.values[region, 2]
      results <- c(pval, is.err, beta.values.estimate, beta.values.estimate - quantile.alpha * beta.values.sd, beta.values.estimate + quantile.alpha * beta.values.sd)
    }else{ results <- c(pval, is.err)}
    
    #Cleaning temporary objects
    rm(score) ; rm(data) ; rm(fit) ; gc()
    return(results)

}
  

  
