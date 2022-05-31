##rbm with GRR
rbm.GRR <- function(genes.maf = Kryukov, size, prev, replicates, 
                    GRR.matrix.del, GRR.matrix.pro = NULL, p.causal = 0.5, p.protect = 0, 
                    same.variant=FALSE, 
                    genetic.model=c("general", "multiplicative", "dominant", "recessive"), select.gene,
                    selected.controls = T, max.maf.causal = 0.01) {
  
  if (nlevels(genes.maf$gene) > 1){ 
    if(missing(select.gene)){
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(genes.maf$gene)[[1]]
    }
    pop.maf <- subset(genes.maf, genes.maf$gene %in% select.gene & genes.maf$maf > 0)$maf
    if(any(subset(genes.maf, genes.maf$gene %in% select.gene)$maf == 0)){
      warning("Some variants have a maf equal to 0 and won't be kept")
      genes.maf <- subset(genes.maf, genes.maf$maf > 0 )
    }
  }else{
    if(any(genes.maf$maf == 0)){
      warning("Some variants have a maf equal to 0 and won't be kept")
      genes.maf <- subset(genes.maf, genes.maf$maf > 0 )
    }
    pop.maf <- genes.maf$maf
  }
  
  ##Check GRR
  if (!is.list(GRR.matrix.del)) {
    if (is.matrix(GRR.matrix.del)) {
      GRR.matrix.del <- list(GRR.matrix.del)
    }else{
      stop("GRR.matrix.del should be a list or a matrix")
    }
  GRR.het <- GRR.matrix.del[[1]]
  }

  genetic.model <- match.arg(genetic.model)

  if (length(GRR.matrix.del) == 1) {
    if (genetic.model == "general") {
      stop("Needs two GRR matrices in the general model")
    }else{
      GRR.homo.alt <- NULL
    }
  }else{
    if (genetic.model == "general") {
      GRR.homo.alt <- GRR.matrix.del[[2]]
    }else{
      warning("Only one GRR matrix needed for this model, only the first one is used")
      GRR.homo.alt <- NULL
    }
  }
  
  ##Same for protective
  if (!is.null(GRR.matrix.pro)) {
    if (!is.list(GRR.matrix.pro)) {
      if (is.matrix(GRR.matrix.pro)) {
        GRR.matrix.pro <- list(GRR.matrix.pro)
      }else{
        stop("GRR.matrix.pro should be a list or a matrix")
      }
    }
    GRR.het.pro <- GRR.matrix.pro[[1]]
    if (length(GRR.matrix.pro) == 1) {
      if (genetic.model == "general") {
        stop("Needs two GRR matrices in the general model")
      }else{
        GRR.homo.alt.pro <- NULL
      }
    }else{
      if (genetic.model == "general") {
        GRR.homo.alt.pro <- GRR.matrix.pro[[2]]
      }else{
        warning("Only one GRR matrix needed for this model, only the first one is used")
        GRR.homo.alt.pro <- NULL
      }
    }
  }else{
    GRR.het.pro <- 1/GRR.het
    GRR.homo.alt.pro <- 1/GRR.homo.alt
  }
  
  ##Check on GRR values
  if (any(GRR.het < 1) | any(GRR.homo.alt < 1)) stop("Matrix of deleterious GRR has GRR values lower than 1")
  if (!is.null(GRR.matrix.pro)) {
    if (any(GRR.het.pro > 1) | any(GRR.homo.alt.pro > 1)) stop("Matrix of protective GRR has GRR values greater than 1")
  }

  ##Order arguments
  GRR.pars <- list(OR.del = GRR.het, OR.pro = GRR.het.pro, p.causal = p.causal, prob.pro = p.protect)
  if (!is.null(GRR.homo.alt)) {
    GRR.homo.alt.pars <- list(OR.del = GRR.homo.alt, OR.pro = GRR.homo.alt.pro)
  }else{
    GRR.homo.alt.pars <- NULL
  }

  ##Check dimensions
  if(!is.null(GRR.homo.alt.pars)){
    if(nrow(GRR.het) != nrow(GRR.homo.alt.pars$OR.del) |  ncol(GRR.het) != ncol(GRR.homo.alt.pars$OR.del)) stop("GRR.het and GRR.homo.alt have different dimensions")
  }
  
  ##Choose the OR function
  if(same.variant){
    variant.function <- OR.matrix.same.variant
  }else{
    variant.function <- OR.matrix
  }
  
  GRR.pars$n.variants <- length(pop.maf)
  nb_snps <- GRR.pars$n.variants * replicates
  nb_inds <- sum(size)
  GRR.pars$maf <- pop.maf
  GRR.pars$maf.threshold <-  max.maf.causal
  x <- new.bed.matrix(nb_inds, nb_snps);
  x@snps$Causal <- ""
  for(b in 1:replicates) {
    GRR.causal <- do.call( variant.function, GRR.pars)
    GRR.het <- GRR.causal$OR
    if(!is.null(GRR.homo.alt.pars)){
      #Give GRR>1 for the same variants as GRR.het
      for(i in 1:length(prev)){
        GRR.homo.alt[i, which(GRR.het[i,]==1)] <- 1
        GRR.homo.alt[i, which(GRR.het[i,]<1)] <- GRR.homo.alt.pars$OR.pro[i, which(GRR.het[i,]<1)]
      }
    }else{
      GRR.homo.alt=NULL
    }    
    MAFS <- genotypic.freq(genes.maf=genes.maf, GRR.het=GRR.het, GRR.homo.alt=GRR.homo.alt, prev=prev, select.gene=select.gene, genetic.model=genetic.model, selected.controls = selected.controls)
  #Check if problems with model
    if(any(MAFS$freq.homo.ref[1,]>1 | MAFS$freq.het[1,]<0 | MAFS$freq.homo.alt[1,]<0)) stop("Impossible genetic model, please change your parametrization")
    .Call("oz_random_filling_bed_matrix_noHW", PACKAGE = "Ravages", x@bed, MAFS$freq.homo.ref, MAFS$freq.het, size, (b-1)*GRR.pars$n.variants)
    #Add flag for causal variants
    if(same.variant) x@snps$Causal[(b-1)*GRR.pars$n.variants + GRR.causal$causal] <- paste("Gpe", 1:length(prev)+1, collapse = "", sep = "")
    else{
      for(i in 1:length(GRR.causal$causal)){
        x@snps$Causal[(b-1)*GRR.pars$n.variants + GRR.causal$causal[[i]]] <- paste0(x@snps$Causal[(b-1)*GRR.pars$n.variants + GRR.causal$causal[[i]]], "Gpe", i+1)
      }
    }
  }
  x@ped$pheno <- factor(rep.int( 1:length(size) - 1, size))
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = GRR.pars$n.variants) )
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x@snps$Causal[which(x@snps$Causal=="")] <- "NonCausal"
  x <- set.stats(x, verbose = FALSE)
  x
}

