
burden.mlogit <- function(x, group = x@ped$pheno, 
                          genomic.region = x@snps$genomic.region, 
                          burden, maf.threshold = 0.01, 
                          ref.level, formula = NULL, data = NULL, get.OR.value=FALSE, alpha=0.05){

  if(!(ref.level %in% levels(group))) 
    stop("'ref.level' is not a level of 'group'")

  if(is.numeric(burden)) {
    if(!is.matrix(burden)){
      stop("Score is not a matrix")
    } else {
      if(ncol(burden) != nlevels(x@snps$genomic.region) | nrow(burden) != nrow(x@ped))
        stop("Score has wrong dimensions")
    }
    score <- burden
  } else if(burden == "CAST"){
    score <- CAST(x, genomic.region, maf.threshold)
  } else if(burden == "WSS"){
    score <- WSS(x, genomic.region)
  } else {
    stop("'burden' should be \"CAST\", \"WSS\", or a matrix of pre-computed burdens");
  }
  score <- as.data.frame(score)

  if(!is.null(data)){
    if(nrow(data) != length(group)){
      stop("'data' has wrong dimensions")
    }
  }
  
  group <- as.factor(group)
  genomic.region <- as.factor(genomic.region)
  
  alt.levels <- levels(group)[levels(group) != ref.level]

  # preparation data / formula
  if(is.null(data)) {
    data <- cbind(ind.pheno = group, score)
  } else {
    data <- as.data.frame(data)
    if(is.null(formula)) 
      formula <- as.formula( paste("~", paste(colnames(data), collapse = "+")))
    data <- cbind(ind.pheno = group, score, data)
  }
  rownames(data) <- NULL # is this useful ??
  data <- mlogit.data(data, varying=NULL, shape="wide", choice="ind.pheno", alt.levels=levels(group))
 
  R <- sapply( levels(genomic.region), 
               function(reg) 
                  run.mlogit(pheno = group, score = score, region = reg, 
                  ref.level = ref.level, alt.levels = alt.levels, formula = formula, data = data, 
                  alpha = alpha, get.OR.value = get.OR.value))

  R <- as.data.frame( t(R) );

  if(get.OR.value)
    colnames(R) <- c("p.value", "is.err", paste("OR", alt.levels, sep="."), paste("l.lower", alt.levels, sep="."), paste("l.upper", alt.levels, sep="."))
  else
    colnames(R) <- c("p.value", "is.err")

  rownames(R) <- levels(genomic.region)
  R
}


run.mlogit <- function(pheno, score, region, ref.level, alt.levels, formula, data, alpha, get.OR.value){
  # Formula for the current region
  if(is.null(formula)) { 
    my.formula <- mFormula(as.formula(paste("ind.pheno ~ 0 |", region)))
  } else {
    z <- as.character(formula)
    if(z[1] != "~" | length(z) != 2) 
      stop("'formula' should be a formula of the form \"~ var1 + var2\"")
    z <- z[2]
    my.formula <- mFormula( as.formula( paste("ind.pheno ~ 0 |", region, " + ", z) ) )
    my.formula.H0 <- mFormula( as.formula( paste("ind.pheno ~ 0 | ", z ) ) )
  }

 
  # Catch errors
  fit <- tryCatch(mlogit(my.formula, data = data, reflevel = ref.level), error = identity, warning = identity)

  if(is(fit, "error")) {
    pval <- NA ; 
    is.err <- 1 ; 
    OR.values <- data.frame("Estimate" = rep(NA, nlevels(pheno)-1), 
                            "sd" = rep(NA, nlevels(pheno)-1)) ; 
    rownames(OR.values) <- paste(alt.levels, "region", sep=":")
  } else {
    if(is.null(formula)) {
      my.model <- summary(fit)
      pval <- as.numeric(my.model$lratio$p.value)
      if(get.OR.value == TRUE)
        OR.values <- my.model$CoefTable
    } else {
      my.model.H0 <- summary(mlogit(my.formula.H0, data = data, reflevel = ref.level))
      my.model.H1 <- summary(fit)
      pval <- pchisq(-2*as.numeric(my.model.H0$logLik) + 2*as.numeric(my.model.H1$logLik), nlevels(pheno)-1, lower.tail=FALSE)
      if(get.OR.value == TRUE) 
        OR.values <- my.model.H1$CoefTable
    }
    is.err <- 0
  }

  quantile.alpha <- qnorm(alpha/2, lower.tail=FALSE)

  if(get.OR.value) 
    results <- c(pval, is.err, 
        as.numeric(exp(OR.values[paste(alt.levels, region, sep=":"),1])), 
        as.numeric(exp(OR.values[paste(alt.levels, region, sep=":"),1]-quantile.alpha*OR.values[paste(alt.levels, region, sep=":"),2])), 
        as.numeric(exp(OR.values[paste(alt.levels, region, sep=":"),1]+quantile.alpha*OR.values[paste(alt.levels, region, sep=":"),2]))
      )
   else
    results <- c(pval, is.err)
  
  return(results)
}
  
