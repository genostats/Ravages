#We only accept matrices of good dimensions
group.mafs.GRR <- function(file.pop.maf = Ravages::Kryukov, GRR, GRR.2=NULL, baseline, model=c("general", "multiplicative", "dominant", "recessive"), select.gene=NULL) {
  #Selection of maf
  if (nlevels(file.pop.maf$gene) > 1) {
    if(is.null(select.gene)) warning("More than one gene in the file")
    pop.maf <- subset(file.pop.maf, file.pop.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- file.pop.maf$maf
  }

  #test dimensions of GRR
  if(nrow(GRR) != length(baseline) | ncol(GRR) != length(pop.maf)) 
    stop("GRR dimensions mismatch")
    
  #Test on GRR.2 only for general model , if not general model: only first GRR matrix used
  if(model=="general"){
    if(is.null(GRR.2)){
      stop("general model needs two GRR values per group")
    }else{
      if(nrow(GRR)!=nrow(GRR.2) | ncol(GRR)!=ncol(GRR.2)){
        stop("GRR and GRR.2 have different dimensions")
      }
    }
  }else{
    if(!is.null(GRR.2)){
      warning("Only one GRR matrix needed for this model, only the first one is used")
    }
  }

  
      
  #Creates matrix for each model
  if(model=="multiplicative"){
    GRR.2 <- GRR**2
  }
  
  if(model=="dominant"){
    GRR.2 <- GRR
  }
  
  #GRR=1 for heterozygous if model=recessive
  if(model=="recessive"){
    GRR.2 <- GRR
    GRR <- matrix(rep(1, ncol(GRR)*nrow(GRR)), nrow=nrow(GRR)) 
  }    

  if(any(baseline < 0) | sum(baseline) > 1)
    stop("Unappropriate baseline values")

  #Frequencies calculation
  p.c <- matrix(rep(0,ncol(GRR)*nrow(GRR)), nrow=nrow(GRR))
  p.t <- numeric(ncol(GRR))
  for (i in 1:ncol(GRR)){
    p.c[,i] <- p.case(pop.maf[i], GRR[,i], GRR.2[,i])
    p.t[i] <- p.tem.GRR(pop.maf[i], GRR[,i], GRR.2[,i], baseline=baseline)
  }
  
  F <- rbind(p.t, p.c)
  if(nrow(F) == 2) 
    rownames(F) <- c("controls","cases")
  else 
    rownames(F) <- c("controls", sprintf("cases_%d", 1:(nrow(F)-1)))
  F
} 
  
  
p.case <- function(p, GRR, GRR.2){
  freq.homo.alt <- (GRR.2*p**2)/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)
  freq.het <- (GRR*2*p*(1-p))/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)
  return(freq.homo.alt + 0.5*freq.het)
}

p.tem.GRR <- function(p, GRR, GRR.2, baseline){
  f <- baseline / (GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2) #Frequence des cas chez les homo de ref
  
  freq.homo.alt <- (p**2 * (1-sum(f*GRR.2))) / (1-sum(baseline))
  freq.het <- (2*p*(1-p) * (1-sum(f*GRR))) / (1-sum(baseline))
  return(freq.homo.alt + 0.5*freq.het)
}



