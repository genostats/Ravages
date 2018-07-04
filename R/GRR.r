group.mafs.GRR <- function(pop.maf, GRR.1, GRR.2=NULL, baseline, model=c("general", "additive", "dominant", "recessive")) {
    if(model=="recessive") GRR.1 <- GRR.2
    if(!is.matrix(GRR.1)) {
    if(length(GRR.1) == length(pop.maf) & length(baseline) == 1)
      GRR.1 <- matrix(GRR.1, nrow = 1)
    else if(length(GRR.1) == length(baseline))
      GRR.1 <- matrix( rep_len(GRR.1, length(pop.maf)*length(GRR.1)), nrow = length(baseline))
    else
      stop("GRR.1 dimension mismatch")
  }

  if(nrow(GRR.1) == length(baseline) & ncol(GRR.1) < length(pop.maf)) 
    GRR <- matrix( rep_len(GRR.1, length(pop.maf)*nrow(GRR.1)), nrow = nrow(GRR.1))
  
  if(nrow(GRR.1) != length(baseline) | ncol(GRR.1) != length(pop.maf)) 
    stop("GRR.1 dimensions mismatch")

  if(length(baseline) != nrow(GRR.1)) 
    baseline <- rep_len(baseline, nrow(GRR.1))
    
  #Test on GRR.2 only for general model  
  if(model=="general"){
    if(!is.null(GRR.2)){
    if(!is.matrix(GRR.2)) {
    if(length(GRR.2) == length(pop.maf) & length(baseline = 1))
      GRR.2 <- matrix(GRR.2, nrow = 1)
    else if(length(GRR.2) == length(baseline))
      GRR.2 <- matrix( rep_len(GRR.2, length(pop.maf)*nrow(GRR.2)), nrow = length(baseline))
    else
      stop("GRR.2 dimension mismatch")
    }
  
  if(nrow(GRR.2) == length(baseline) & ncol(GRR.2) < length(pop.maf)) 
    GRR <- matrix( rep_len(GRR.2, length(pop.maf)*nrow(GRR.2)), nrow = nrow(GRR.2))
  
  if(nrow(GRR.1) != nrow(GRR.2) | ncol(GRR.1) != ncol(GRR.2))
    stop("GRR.1 and GRR.2 have different dimensions")
  }
  else{
    stop("general model needs two GRR value per group")
  }}
  
  #Creates matrix for each model
  if(model=="additive"){
    if(!is.null(GRR.2)) stop("only one GRR value per group needed")
    GRR.2 <- GRR.1**2
  }
  
  if(model=="dominant"){
    if(!is.null(GRR.2)) stop("only one GRR value per group needed")
    GRR.2 <- GRR.1
  }
  
  if(model=="recessive"){
    GRR.1 <- matrix(rep(1, nrow(GRR.2)*ncol(GRR.2)), nrow=nrow(GRR.2))
  }  

  if(any(baseline < 0) | sum(baseline) > 1)
    stop("Unappropriate baseline values")

  #Frequencies calculation
  p.c <- matrix(rep(0,ncol(GRR.1)*nrow(GRR.1)), nrow=nrow(GRR.1))
  p.t <- numeric(ncol(GRR.1))
  for (i in 1:ncol(GRR.1)){
    p.c[,i] <- p.case(pop.maf[i], GRR.1[,i], GRR.2[,i])
    p.t[i] <- p.tem.GRR(pop.maf[i], GRR.1[,i], GRR.2[,i], baseline=baseline)
  }
  
  F <- rbind(p.t, p.c)
  if(nrow(F) == 2) 
    rownames(F) <- c("controls","cases")
  else 
    rownames(F) <- c("controls", sprintf("cases_%d", 1:(nrow(F)-1)))
  F
} 
  
  
p.case <- function(p, GRR.1, GRR.2){
  freq.homo.alt <- (GRR.2*p**2)/(GRR.2*p**2 + GRR.1*2*p*(1-p) + (1-p)**2)
  freq.het <- (GRR.1*2*p*(1-p))/(GRR.2*p**2 + GRR.1*2*p*(1-p) + (1-p)**2)
  return(freq.homo.alt + 0.5*freq.het)
}

p.tem.GRR <- function(p, GRR.1, GRR.2, baseline){
  f <- baseline / (GRR.2*p**2 + GRR.1*2*p*(1-p) + (1-p)**2) #Frequence des cas chez les homo de ref
  
  freq.homo.alt <- (p**2 * (1-sum(f*GRR.2))) / (1-sum(baseline))
  freq.het <- (2*p*(1-p) * (1-sum(f*GRR.1))) / (1-sum(baseline))
  return(freq.homo.alt + 0.5*freq.het)
}
