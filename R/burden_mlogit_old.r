burden.mlogit.old <- function(x, group = x@ped$pheno,
                          genomic.region = x@snps$genomic.region, 
                          burden, maf.threshold = 0.01, 
                          ref.level, formula, data, get.OR.value=FALSE, alpha=0.05){
  
  if(!is.factor(group)) group <- as.factor(group)

  if(is.numeric(ref.level)) ref.level <- as.character(ref.level)

  if(!(ref.level %in% levels(group))) 
    stop("'ref.level' is not a level of 'group'")

  if(is.numeric(burden)) {
    if(!is.matrix(burden)){
      stop("Score is not a matrix")
    } else {
      if(is.null(colnames(burden))){ 
        colnames(burden) <- make.names(1:ncol(burden))
      }
    }
    score <- burden
  } else{
    if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
    genomic.region <- droplevels(genomic.region)
    if(missing(x)) stop("a bed.matrix 'x' is needed to compute the score")

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

  if(missing(data)) 
    data <- NULL
  if(missing(formula))
    formula <- NULL

  if(!is.null(data)) {
    if(nrow(data) != length(group)){
      stop("'data' has wrong dimensions")
    }
  }

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
  data.reg <- dfidx(data.reg, varying=NULL, shape="wide", choice="ind.pheno")
 
  R <- sapply( names(score), 
               function(reg) 
                  run.mlogit(pheno = group, score = score, region = reg, 
                  ref.level = ref.level, alt.levels = alt.levels, formula = formula, data = data.reg, 
                  alpha = alpha, get.OR.value = get.OR.value))

  R <- as.data.frame( t(R) );

  if(get.OR.value)
    colnames(R) <- c("p.value", "is.err", paste("OR", alt.levels, sep="."), paste("l.lower", alt.levels, sep="."), paste("l.upper", alt.levels, sep="."))
  else
    colnames(R) <- c("p.value", "is.err")

  rownames(R) <- old.names

  return(R)
}


run.mlogit <- function(pheno, score, region, ref.level, alt.levels, formula, data, alpha, get.OR.value){
  # Formula for the current region
  if(is.null(formula)) { 
    my.formula <- Formula(as.formula(paste("ind.pheno ~ 0 |", region)))
  } else {
    z <- as.character(formula)
    if(z[1] != "~" | length(z) != 2) 
      stop("'formula' should be a formula of the form \"~ var1 + var2\"")
    z <- z[2]
    my.formula <- Formula( as.formula( paste("ind.pheno ~ 0 |", region, " + ", z) ) )
    my.formula.H0 <- Formula( as.formula( paste("ind.pheno ~ 0 | ", z ) ) )
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
                 as.numeric(exp(OR.values[paste(region, alt.levels, sep=":"),1])), 
                 as.numeric(exp(OR.values[paste(region, alt.levels, sep=":"),1]-quantile.alpha*OR.values[paste(region, alt.levels, sep=":"),2])), 
                 as.numeric(exp(OR.values[paste(region, alt.levels, sep=":"),1]+quantile.alpha*OR.values[paste(region, alt.levels, sep=":"),2]))
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
