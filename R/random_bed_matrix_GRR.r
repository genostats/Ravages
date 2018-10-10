##random.bed.matrix with GRR
random.bed.matrix.GRR <- function(file.pop.maf = Kryukov, size, baseline, replicates, GRR.matrix, GRR.matrix.pro=NULL, prop.del = 0.5, prop.pro = 0, 
								  same.variant=c(FALSE, TRUE), fixed.variant.prop = c(TRUE, FALSE), 
								  genetic.model=c("general", "multiplicative", "dominant", "recessive"), select.gene=NULL) {
  
  if (nlevels(file.pop.maf$gene) > 1){ 
    if(is.null(select.gene)){
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(file.pop.maf$gene)[[1]]
    }
    pop.maf <- subset(file.pop.maf, file.pop.maf$gene %in% select.gene)$maf
  }else{
    pop.maf <- file.pop.maf$maf
  }
  
  ##Check GRR
  if (!is.list(GRR.matrix)) {
    if (is.matrix(GRR.matrix)) {
      GRR.matrix <- list(GRR.matrix)
    }else{
      stop("GRR.matrix should be a list or a matrix")
    }
  GRR <- GRR.matrix[[1]]
  }
  if (length(GRR.matrix) == 1) {
    if (genetic.model == "general") {
      stop("Needs two GRR matrices in the general model")
    }else{
      GRR.2 <- NULL
    }
  }else{
    if (genetic.model == "general") {
      GRR.2 <- GRR.matrix[[2]]
    }else{
      warning("Only one GRR matrix needed for this model, only the first one is used")
      GRR.2 <- NULL
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
    GRR.pro <- GRR.matrix.pro[[1]]
    if (length(GRR.matrix.pro) == 1) {
      if (genetic.model == "general") {
        stop("Needs two GRR matrices in the general model")
      }else{
        GRR.2.pro <- NULL
      }
    }else{
      if (genetic.model == "general") {
        GRR.2.pro <- GRR.matrix.pro[[2]]
      }else{
        warning("Only one GRR matrix needed for this model, only the first one is used")
        GRR.2.pro <- NULL
      }
    }
  }else{
    GRR.pro <- 1/GRR
    GRR.2.pro <- 1/GRR.2
  }
  
  ##Check on GRR values
  if (any(GRR < 1) | any(GRR.2 < 1)) stop("Matrix of deleterious GRR has GRR values lower than 1")
  if (!is.null(GRR.matrix.pro)) {
    if (any(GRR.pro > 1) | any(GRR.2.pro > 1)) stop("Matrix of protective GRR has GRR values greater than 1")
  }

  ##Order arguments
  GRR.pars <- list(OR.del = GRR, OR.pro = GRR.pro, prob.del = prop.del, prob.pro = prop.pro)
  if (!is.null(GRR.2)) {
    GRR.2.pars <- list(OR.del = GRR.2, OR.pro = GRR.2.pro)
  }else{
    GRR.2.pars <- NULL
  }

  ##Check dimensions
  if(!is.null(GRR.2.pars)){
    if(nrow(GRR) != nrow(GRR.2.pars$OR.del) |  ncol(GRR) != ncol(GRR.2.pars$OR.del)) stop("GRR and GRR.2 have different dimensions")
  }
  
  ##Choose the OR function
  if(fixed.variant.prop){
    variant.function <- ifelse(same.variant == FALSE, OR.matrix.fix, OR.matrix.fix.same.variant)
  }else{
    variant.function <- ifelse(same.variant == FALSE, OR.matrix, OR.matrix.same.variant)
  }
  
  GRR.pars$n.variants <- length(pop.maf)
  nb_snps <- GRR.pars$n.variants * replicates
  nb_inds <- sum(size)
  x <- new.bed.matrix(nb_inds, nb_snps);
  for(b in 1:replicates) {
    GRR <- do.call( variant.function, GRR.pars)
    if(!is.null(GRR.2.pars)){
      #Give GRR>1 for the same variants as GRR
      for(i in 1:length(baseline)){
        GRR.2[i, which(GRR[i,]==1)] <- 1
        GRR.2[i, which(GRR[i,]<1)] <- GRR.2.pars$OR.pro[i, which(GRR[i,]<1)]
      }
    }else{
      GRR.2=NULL
    }    
    MAFS <- genotypic.freq(file.pop.maf=file.pop.maf, GRR=GRR, GRR.2=GRR.2, baseline=baseline, select.gene=select.gene, genetic.model=genetic.model)
    .Call("oz_random_filling_bed_matrix_noHW", PACKAGE = "Ravages", x@bed, MAFS$freq.homo.ref, MAFS$freq.het, size, (b-1)*GRR.pars$n.variants)
  }
  x@ped$pheno <- rep.int( 1:length(size) - 1, size)
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = GRR.pars$n.variants) )
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x <- set.stats(x, verbose = FALSE)
  x
}

