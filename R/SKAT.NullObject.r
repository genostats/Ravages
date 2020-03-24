SKAT.NullObject <- function(group, data, formula){
  if(!is.factor(group))
    stop("'group' is not a factor")
    
  ref.level <- levels(group)[1]

  if(missing(data)) {
    a <- table(group)/length(group)
    Pi.data <- matrix( a, ncol = nlevels(group), nrow = length(group), byrow = TRUE)
    X <- matrix(1, nrow = length(group), ncol = 1) 
  }else{
    if(!is.matrix(data)) stop("'data' should be a matrix")
    if(is.null(colnames(data))) colnames(data) <- sprintf("C%0*d", log10(ncol(data)) + 1, 1:ncol(data))
    if(missing(formula)) formula <- as.formula(paste("~", paste(colnames(data), collapse = "+")))
    Pi.data <- Pi.matrix(group, data, formula, ref.level)
    X <- cbind(1, data[, all.vars(formula), drop=FALSE]) 
  }
  
  return(list(Pi.data = Pi.data, X = X, group = group))
}
