#We only accept matrices of good dimensions
group.mafs.GRR <- function(file.pop.maf = Ravages::Kryukov, GRR, GRR.2=NULL, baseline, model=c("general", "multiplicative", "dominant", "recessive"), select.gene) {
  #Selection of maf
  if (nlevels(file.pop.maf$gene) > 1) {
    if(missing(select.gene)){
      warning("More than one gene in the file")
      select.gene <- levels(file.pop.maf$gene)[[1]]
    }
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
  homo.ref <- het <- homo.alt <- matrix(rep(0, ncol(GRR)*(nrow(GRR)+1)), nrow=nrow(GRR)+1)
  p.c <- matrix(rep(0,ncol(GRR)*nrow(GRR)), nrow=nrow(GRR))
  p.t <- numeric(ncol(GRR))
  for (i in 1:ncol(GRR)){
    freq.case <- p.case(pop.maf[i], GRR[,i], GRR.2[,i])
    freq.controls <- p.tem.GRR(pop.maf[i], GRR[,i], GRR.2[,i], baseline=baseline)

    homo.ref[,i] <- c(freq.controls["freq.homo.ref",], freq.case["freq.homo.ref",])
    het[,i] <- c(freq.controls["freq.het",], freq.case["freq.het",])
    homo.alt[,i] <- c(freq.controls["freq.homo.alt",], freq.case["freq.homo.alt",])
  }
  
  if(nrow(homo.ref) == 2) 
    rownames(homo.ref) <- rownames(het) <- rownames(homo.alt) <- c("controls","cases")
  else 
    rownames(homo.ref) <- rownames(het) <- rownames(homo.alt) <- c("controls", sprintf("cases_%d", 1:(nrow(homo.ref)-1)))
  return(list(freq.homo.ref = homo.ref, freq.het = het, freq.homo.alt = homo.alt))
} 
  
  
p.case <- function(p, GRR, GRR.2){
  freq.homo.alt <- (GRR.2*p**2)/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)
  freq.het <- (GRR*2*p*(1-p))/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)
  freq.homo.ref <- 1-freq.homo.alt-freq.het
  return(t(data.frame(freq.homo.ref, freq.het, freq.homo.alt)))
}

p.tem.GRR <- function(p, GRR, GRR.2, baseline){
  f <- baseline / (GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2) #Frequence des cas chez les homo de ref
  
  freq.homo.alt <- (p**2 * (1-sum(f*GRR.2))) / (1-sum(baseline))
  freq.het <- (2*p*(1-p) * (1-sum(f*GRR))) / (1-sum(baseline))
  freq.homo.ref <- 1-freq.homo.alt-freq.het
  return(t(data.frame(freq.homo.ref, freq.het, freq.homo.alt)))
}



