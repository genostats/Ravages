burden.continuous <- function(x, pheno = x@ped$pheno, genomic.region = x@snps$genomic.region, burden, maf.threshold = 0.01, formula, data){
 
  if(!missing(x)){
    if(length(genomic.region)==0){
      warning("No 'genomic region' given, all variants will be analysed in the same testing unit")
      genomic.region <- rep("UniqRegion", ncol(x))
    }
  }  
  if(!is.numeric(pheno)) stop("'pheno' should be a numeric vector")

  if(is.numeric(burden)) {
    if(!is.matrix(burden)){
      stop("Score is not a matrix")
    } else {
      if(is.null(colnames(burden))){ 
        colnames(burden) <- make.names(1:ncol(burden))
      }
    }
    score <- burden
  } else if(burden == "CAST"){
    if(missing(x)) stop("a bed.matrix 'x' is needed to compute CAST score")
    score <- CAST(x, genomic.region, maf.threshold)
  } else if(burden == "WSS"){
    if(missing(x)) stop("a bed.matrix 'x' is needed to compute WSS score")
    score <- WSS(x, genomic.region)
  } else {
    stop("'burden' should be \"CAST\", \"WSS\", or a matrix of pre-computed burdens");
  }
  score <- as.data.frame(score)
  # to ensure syntactically correct formulas
  old.names <- colnames(score)
  names(score) <- make.names(names(score))

  if(missing(data)) 
    data <- NULL
  if(missing(formula)){
    formula <- NULL
  }

  if(!is.null(data)) {
    if(is.null(colnames(data))) colnames(data) <- sprintf("C%0*d", log10(ncol(data))+1, 1:ncol(data))
    if(nrow(data) != length(pheno)){
      stop("'data' has wrong dimensions")
    }
  }

  # preparation data / formula
  if(is.null(data)) {
    data.reg <- cbind(ind.pheno = pheno, score)
  } else {
    data.reg <- as.data.frame(data)
    if(is.null(formula))  
      formula <- as.formula( paste("~", paste(colnames(data), collapse = "+")))
    data.reg <- cbind(ind.pheno = pheno, score, data.reg)
  }


  R <- sapply( names(score), function(reg) run.burden.continuous(pheno = pheno, score = score, region = reg, formula = formula, data = data.reg))

  R <- as.data.frame( t(R) );

  colnames(R) <- c("p.value", "is.err")

  rownames(R) <- old.names

  return(R)
}


run.burden.continuous <- function(pheno, score, region, formula, data){
  # Formula for the current region
  if(is.null(formula)) { 
    my.formula <- as.formula(paste("ind.pheno ~ ", region))
  } else {
    z <- as.character(formula)
    if(z[1] != "~" | length(z) != 2) 
      stop("'formula' should be a formula of the form \"~ var1 + var2\"")
    z <- z[2]
    my.formula <- as.formula( paste("ind.pheno ~ ", region, " + ", z) ) 
  }

 
  # Catch errors
  fit <- tryCatch(lm(my.formula, data = data), error = identity, warning = identity)

  if(is(fit, "error")) {
    pval <- NA ; 
    is.err <- 1 ; 
  } else {
    my.model <- summary(fit)
    pval <- my.model$coefficients[region, 4]
    is.err <- 0
  }

  return(c(pval, is.err))
}
  

  