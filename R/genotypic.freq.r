#We only accept matrices of good dimensions
genotypic.freq <- function(genes.maf = Kryukov, GRR, GRR.2=NULL, baseline, genetic.model=c("general", "multiplicative", "dominant", "recessive"), select.gene=NULL) {

  #Test if a good genetic model is given
  if(!(genetic.model %in% c("general", "multiplicative", "dominant", "recessive"))) stop("Wrong genetic.model")

  #Selection of maf
  if (nlevels(genes.maf$gene) > 1) {
    if(is.null(select.gene)){
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(genes.maf$gene)[[1]]
    }
    pop.maf <- subset(genes.maf, genes.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- genes.maf$maf
  }

  #test dimensions of GRR
  if(nrow(GRR) != length(baseline) | ncol(GRR) != length(pop.maf)) 
    stop("GRR dimensions mismatch")
    
  #Test on GRR.2 only for general model , if not general model: only first GRR matrix used
  if(genetic.model=="general"){
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
  if(genetic.model=="multiplicative"){
    GRR.2 <- GRR**2
  }
  
  if(genetic.model=="dominant"){
    GRR.2 <- GRR
  }
  
  #GRR=1 for heterozygous if genetic.model=recessive
  if(genetic.model=="recessive"){
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

    homo.ref[,i] <- c(freq.controls[3], freq.case[,3])
    het[,i] <- c(freq.controls[2], freq.case[,2])
    homo.alt[,i] <- c(freq.controls[1], freq.case[,1])
  }
  
  if(nrow(homo.ref) == 2) 
    rownames(homo.ref) <- rownames(het) <- rownames(homo.alt) <- c("controls","cases")
  else 
    rownames(homo.ref) <- rownames(het) <- rownames(homo.alt) <- c("controls", sprintf("cases_%d", 1:(nrow(homo.ref)-1)))
  return(list(freq.homo.ref = homo.ref, freq.het = het, freq.homo.alt = homo.alt))
} 
  
  
p.case <- function(p, GRR, GRR.2){
  r <- matrix( 0.0, ncol = 3, nrow = length(GRR))
  r[,1] <- (GRR.2*p**2)/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)     # freq homo alt
  r[,2] <- (GRR*2*p*(1-p))/(GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2)  # freq het
  r[,3] <- 1-r[,1]-r[,2]                                            # freq homo ref
  return(r)
}

p.tem.GRR <- function(p, GRR, GRR.2, baseline){
  f <- baseline / (GRR.2*p**2 + GRR*2*p*(1-p) + (1-p)**2) #Frequence des cas chez les homo de ref
  r <- numeric(3)
  r[1] <- (p**2 * (1-sum(f*GRR.2))) / (1-sum(baseline))           # freq homo alt
  r[2] <- (2*p*(1-p) * (1-sum(f*GRR))) / (1-sum(baseline))        # freq het
  r[3] <- 1-r[1]-r[2]                                             # freq homo ref
  return(r)
}



