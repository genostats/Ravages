compute.GRR.matrix <- function(file.pop.maf=Ravages::Kryukov, n.case.groups=2, GRR=c("constant", "SKAT", "variable"), GRR.value=NULL, GRR.formula=NULL, GRR.multiplicative.factor=2, select.gene="R1"){
  ##Select MAF from the file given by the user  
  if(nlevels(file.pop.maf$gene)>1){
    if(is.null(select.gene)) warning("More than one gene in the file")
    pop.maf <- subset(file.pop.maf, file.pop.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- file.pop.maf$maf
  }
                  
  n.variants <- length(pop.maf)
                      
  ##Same GRR for all variants
  if(GRR=="constant"){
    if(is.null(GRR.value)) stop("Needs GRR value for variants")
    if(length(GRR.value)>1){
      warning("Only one GRR value needed, only the first value will be used")
      GRR.value <- GRR.value[1]
    }
    if(is.null(GRR.multiplicative.factor)) stop("Needs multiplicative factors")
    #Need n.groups-1 multiplicative factor: multiplication between group of case 1 and other groups of cases
    if(length(GRR.multiplicative.factor)!=(n.case.groups-1)) stop("Wrong number of multiplicative factors")
    GRR.matrix <- matrix(rbind(rep(GRR.value, n.variants), t(sapply(GRR.multiplicative.factor, function(z) z*rep(GRR.value, n.variants)))), nrow=n.case.groups)
  }
                                                                                                  
  if(GRR=="SKAT"){
    if(is.null(pop.maf)) stop("Needs MAF in the population to compute GRR")
    if(is.null(GRR.multiplicative.factor)) stop("Needs multiplicative factors")
    if(length(GRR.multiplicative.factor)!=(n.case.groups-1)) stop("Wrong number of multiplicative factors")
    GRR.matrix <- matrix(rbind(exp(0.402*abs(log10(pop.maf))), t(sapply(GRR.multiplicative.factor, function(z) z*exp(0.402*abs(log10(pop.maf)))))), nrow=n.case.groups)
  }
                                                                              
  if(GRR=="variable"){
    if(is.null(pop.maf)) stop("Needs MAF in the population to compute GRR")
    if(is.null(GRR.formula)) stop("Needs a formula to compute GRR")
    if(is.null(GRR.multiplicative.factor)) stop("Needs multiplicative factors")
    if(length(GRR.multiplicative.factor)!=(n.case.groups-1)) stop("Wrong number of multiplicative factors")
    GRR.matrix <- matrix(rbind(do.call(GRR.formula, list(pop.maf)), t(sapply(GRR.multiplicative.factor, function(z) z*do.call(GRR.formula, list(pop.maf))))), nrow=n.case.groups, byrow=TRUE)
  }
  return(GRR.matrix)
}
