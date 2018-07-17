#We only accept matrices of good dimensions
group.mafs.GRR <- function(pop.maf, GRR, GRR.2=NULL, baseline, model=c("general", "additive", "dominant", "recessive")) {
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
  if(model=="additive"){
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




##random.bed.matrix with GRR
random.bed.matrix.GRR <- function(pop.maf, size, baseline, replicates, GRR.pars, GRR.pars.2=NULL, OR.function, model=c("general", "additive", "dominant", "recessive")) {
  GRR.pars$n.variants <- length(pop.maf)
  nb_snps <- GRR.pars$n.variants * replicates
  nb_inds <- sum(size)
  x <- oz:::new.bed.matrix(nb_inds, nb_snps);
  for(b in 1:replicates) {
    GRR <- do.call( OR.function, GRR.pars)
    if(!is.null(GRR.pars.2)){
    #Test if good dimensions
      if(nrow(GRR) != nrow(GRR.pars.2$OR.del) |  ncol(GRR) != ncol(GRR.pars.2$OR.del)) stop("GRR and GRR.2 have different dimensions")
    #Protective OR: 1/OR.del if no argument for OR.pro
      if(is.null(GRR.pars.2$OR.pro)) GRR.pars.2$OR.pro <- 1/GRR.pars.2$OR.del
      GRR.2 <- GRR.pars.2$OR.del
      #Give GRR>1 for the same variants as GRR
      for(i in 1:length(baseline)){
        GRR.2[i, which(GRR[i,]==1)] <- 1
        GRR.2[i, which(GRR[i,]<1)] <- GRR.pars.2$OR.pro[i, which(GRR[i,]<1)]
      }
    }else{
      GRR.2=NULL
    }    
    MAFS <- Ravages:::group.mafs.GRR(pop.maf=pop.maf, GRR=GRR, GRR.2=GRR.2, baseline=baseline, model=model)
    .Call("oz_random_filling_bed_matrix", PACKAGE = "oz", x@bed, MAFS, size, (b-1)*GRR.pars$n.variants)
  }
  x@ped$pheno <- rep.int( 1:length(size) - 1, size)
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = GRR.pars$n.variants) )
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x <- set.stats(x, verbose = FALSE)
  x
}

