
burden.mlogit <- function(x, group = x@ped$pheno,
                          genomic.region = x@snps$genomic.region, 
                          burden, maf.threshold = 0.01, 
                          ref.level, formula=NULL, data=NULL, get.OR.value=FALSE, alpha=0.05){

  group <- if(!is.factor(group)) as.factor(group)

  if(is.numeric(ref.level)) ref.level <- as.character(ref.level)

  if(!(ref.level %in% levels(group))) 
    stop("'ref.level' is not a level of 'group'")

  if(is.numeric(burden)) {
    if(!is.matrix(burden)){
      stop("Score is not a matrix")
    } else {
      #If burden is a matrix: no need to specify genomic.region  
      if(is.null(colnames(burden))){ 
        genomic.region <- factor(paste("genomic.region", 1:ncol(burden), sep="."))
      }else{
        genomic.region <- as.factor(colnames(burden))
      }
      if(ncol(burden) != nlevels(genomic.region) | nrow(burden) != length(group))
        stop("Score has wrong dimensions")
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

  if(!is.null(data)){
    if(nrow(data) != length(group)){
      stop("'data' has wrong dimensions")
    }
  }

  genomic.region <- as.factor(genomic.region)
  
  alt.levels <- levels(group)[levels(group) != ref.level]

  # preparation data / formula
  if(is.null(data)) {
    data.reg <- cbind(ind.pheno = group, score)
  } else {
    data.reg <- as.data.frame(data)
    if(is.null(formula)) 
      formula <- as.formula( paste("~", paste(colnames(data), collapse = "+")))
    data.reg <- cbind(ind.pheno = group, score, data.reg)
  }
  rownames(data.reg) <- NULL # is this useful ??
  data.reg <- mlogit.data(data.reg, varying=NULL, shape="wide", choice="ind.pheno", alt.levels=levels(group))
 
  R <- sapply( levels(genomic.region), 
               function(reg) 
                  run.mlogit(pheno = group, score = score, region = reg, 
                  ref.level = ref.level, alt.levels = alt.levels, formula = formula, data = data.reg, 
                  alpha = alpha, get.OR.value = get.OR.value))

  R <- as.data.frame( t(R) );

  if(get.OR.value)
    colnames(R) <- c("p.value", "is.err", paste("OR", alt.levels, sep="."), paste("l.lower", alt.levels, sep="."), paste("l.upper", alt.levels, sep="."))
  else
    colnames(R) <- c("p.value", "is.err")

  rownames(R) <- levels(genomic.region)

  return(R)
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
  


run.mlogit.two.models <- function(large.group, combined.group, score, region, ref.level, alt.levels.large, alt.levels.combined, formula, data.group, data.combined){
  #Check if same reference level
  if(sum(large.group==ref.level) != sum(combined.group==ref.level)) stop("Not the same reference group in the two models")
  
  # Formula for the current region
  if(is.null(formula)) {
    my.formula.combined.gpe <- mFormula(as.formula(paste("combined.group ~ 0 |", region)))
    my.formula.large.gpe <- mFormula(as.formula(paste("large.group ~ 0 |", region)))
  } else {
    z <- as.character(formula)
    if(z[1] != "~" | length(z) != 2)
      stop("'formula' should be a formula of the form \"~ var1 + var2\"")
    z <- z[2]
    my.formula.combined.gpe <- mFormula( as.formula( paste("combined.group ~ 0 |", region, " + ", z) ) )
    my.formula.large.gpe <- mFormula( as.formula( paste("large.group ~ 0 | ", z ) ) )
  }

  #Catch errors
  fit.large.gpe <- tryCatch(mlogit(my.formula.large.gpe, data = data.group, reflevel = ref.level), error = identity, warning = identity)
  fit.combined.gpe <- tryCatch(mlogit(my.formula.combined.gpe, data = data.combined, reflevel = ref.level), error = identity, warning = identity)
  
  if(is(fit.large.gpe, "error") | is(fit.combined.gpe, "error")) {
    pval <- NA ;
    is.err <- 1 ;
  }else {
    my.model.combined.gpe <- summary(fit.combined.gpe)
    my.model.large.gpe <- summary(fit.large.gpe)
    pval <- pchisq(-2*as.numeric(my.model.combined.gpe$logLik) + 2*as.numeric(my.model.large.gpe$logLik), nlevels(large.group)-nlevels(combined.group), lower.tail=FALSE)
    is.err <- 0
  }
  
  results <- c(pval, is.err)
  return(results)
}
