rbm.GRR.power <- function(genes.maf = Kryukov, size = c(500, 500), prev = 0.01, GRR.matrix.del, GRR.matrix.pro = NULL, p.causal = 0.5, p.protect = 0, same.variant = FALSE, genetic.model = c("multiplicative", "general", "dominant", "recessive"), select.gene, alpha = 2.5e-6, selected.controls = TRUE, power.type = c("simulations", "theoretical"), verbose = TRUE, RVAT = c("CAST", "WSS", "SKAT"), SKAT.method = c("permutations", "theoretical"), max.maf.causal = 0.01, maf.filter = max.maf.causal, replicates = 1000, cores = 10){
  power.type <- match.arg(power.type)
  SKAT.method <- match.arg(SKAT.method)
  genetic.model <- match.arg(genetic.model)
  if(power.type != "simulations" & power.type != "theoretical") stop("'power.type' should be 'simulations' or 'theoretical'")
  if(power.type=="simulations" & !("CAST" %in% RVAT | "WSS" %in% RVAT | "SKAT" %in% RVAT)) stop("Power calculations are only available for 'CAST', 'WSS' and 'SKAT'")
  
  if(missing(select.gene) & nlevels(genes.maf$gene)>1){
    warning("More than one gene in the file, only the first one is used")
    select.gene <- levels(genes.maf$gene)[[1]]
  }
  
  #*With simulations*
  if(power.type=="simulations"){
    #Simulations using rbm.GRR
    x <- rbm.GRR(genes.maf = genes.maf, size = size, prev = prev, GRR.matrix.del = GRR.matrix.del, GRR.matrix.pro = GRR.matrix.pro, p.causal = p.causal, p.protect = p.protect, select.gene = select.gene, same.variant = same.variant, genetic.model = genetic.model, replicates = replicates, selected.controls = selected.controls, max.maf.causal = max.maf.causal)
    pow.names <- c()
    
    #Filtering to keep only rare variants
    x <- filter.rare.variants(x, maf.threshold = maf.filter, min.nb.snps = 2)
    
    ##RVAT
    if("CAST" %in% RVAT | "WSS" %in% RVAT){
      H0.burden <- NullObject.parameters(x@ped$pheno, RVAT = "burden", pheno.type = "cat", ref.level = 0)
      if("CAST" %in% RVAT){
        x.CAST <- burden(x, H0.burden, burden = "CAST", verbose = verbose, maf.threshold = 0.5, cores = cores)
        x.CAST.pow <- sapply(alpha, function(z) mean(x.CAST$p.value<z, na.rm = T))
        pow.names <- c(pow.names, "CAST")
      }else{
        x.CAST.pow <- NULL
      }
      if("WSS" %in% RVAT){
        x.WSS <- burden(x, H0.burden, burden = "WSS", verbose = verbose, cores = cores)
        x.WSS.pow <- sapply(alpha, function(z) mean(x.WSS$p.value<z, na.rm = T))
        pow.names <- c(pow.names, "WSS")
      }else{
        x.WSS.pow <- NULL
      }
    }else{
      x.CAST.pow <- NULL
      x.WSS.pow <- NULL
    }
    if("SKAT" %in% RVAT){
      H0.SKAT <- NullObject.parameters(x@ped$pheno, RVAT = "SKAT", pheno.type = "cat")
      x.SKAT <- SKAT(x, H0.SKAT, verbose = verbose, cores = cores, get.moments = SKAT.method)
      x.SKAT.pow <- sapply(alpha, function(z) mean(x.SKAT$p.value<z, na.rm = T))
      pow.names <- c(pow.names, "SKAT")
    }else{
      x.SKAT.pow <- NULL
    }
  pow <- rbind(x.CAST.pow, x.WSS.pow, x.SKAT.pow)
  rownames(pow) <- pow.names
  colnames(pow) <- alpha 
  } else{
    pow <- CAST.theoretical(genes.maf = genes.maf, size = size, prev = prev, replicates = replicates, GRR.matrix.del = GRR.matrix.del, GRR.matrix.pro = GRR.matrix.pro, p.causal = p.causal, p.protect = p.protect, same.variant = same.variant, genetic.model = genetic.model, select.gene = select.gene, selected.controls = selected.controls, alpha = alpha, maf.threshold = max.maf.causal)
  }
  pow
}
  
  
  
CAST.theoretical <- function(genes.maf = Kryukov, size, prev, replicates, GRR.matrix.del, GRR.matrix.pro, p.causal, p.protect, same.variant, genetic.model, select.gene, selected.controls, alpha, maf.threshold){
  
  #Checks
  if (nlevels(genes.maf$gene) > 1){ 
    if(missing(select.gene)){
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(genes.maf$gene)[[1]]
    }
    genes.maf <- subset(genes.maf, genes.maf$gene %in% select.gene)
    #Keep only rare variants
    var.rare <- which(genes.maf$maf <= maf.threshold)
    pop.maf <- genes.maf$maf[var.rare]
  }else{
    #Keep only rare variants
    var.rare <- which(genes.maf$maf <= maf.threshold)
    pop.maf <- genes.maf$maf[var.rare]
  }
  
  ##Check GRR
  if (!is.list(GRR.matrix.del)) {
    if (is.matrix(GRR.matrix.del)) {
      GRR.matrix.del <- list(GRR.matrix.del)
    }else{
      stop("GRR.matrix.del should be a list or a matrix")
    }
  }
  #Keep only rare variants
  GRR.matrix.del <- lapply(GRR.matrix.del, function(z) z[,var.rare])
  GRR.het <- GRR.matrix.del[[1]]

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
    #Keep only rare variants
    GRR.matrix.pro <- lapply(GRR.matrix.pro, function(z) z[,var.rare])
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
  nb_inds <- sum(size)
  GRR.pars$maf <- pop.maf
  GRR.pars$maf.threshold <- maf.threshold

  #Compute power
  pow <- matrix(NA, ncol = length(alpha), nrow = replicates)
  for (b in 1:replicates) {
      GRR.het <- do.call( variant.function, GRR.pars)$OR
      if(!is.null(GRR.homo.alt.pars)){
        #Give GRR>1 for the same variants as GRR.het
        for(i in 1:length(prev)){
          GRR.homo.alt[i, which(GRR.het[i,]==1)] <- 1
          GRR.homo.alt[i, which(GRR.het[i,]<1)] <- GRR.homo.alt.pars$OR.pro[i, which(GRR.het[i,]<1)]
        } 
      }else{
        if (genetic.model == "multiplicative") {
          GRR.homo.alt <- GRR.het^2
        }
        if (genetic.model == "dominant") {
          GRR.homo.alt <- GRR.het
        }
        if (genetic.model == "recessive") {
          GRR.homo.alt <- GRR.het
          GRR.het <- matrix(rep(1, ncol(GRR.het) * nrow(GRR.het)), nrow = nrow(GRR.het))
        }
      }
         
      P0SC <- sapply(1:length(size[-1]), function(group) prod((1-pop.maf)**2/(GRR.homo.alt[group,]*pop.maf**2+GRR.het[group,]*2*pop.maf*(1-pop.maf)+(1-pop.maf)**2)))
      P1SC <- 1-P0SC
      #In the controls groups
      if(selected.controls){
        PCo <- 1-sum(prev)
        P0ST <- (prod((1-pop.maf)**2)-sum(P0SC*prev))/PCo
      }else{
        P0ST <- prod((1-pop.maf)**2)
      }
      P1ST <- 1-P0ST
      #In all groups of individuals: controls and cases
      P0SG <- c(P0ST, P0SC) ; P1SG <- c(P1ST, P1SC)
      #Probability of the two scores
      P0 <- sum(P0SG*(size/nb_inds))
      P1 <- 1-P0
      #Non-centrality parameter
      ncp <- nb_inds * (sum(size*(P0SG-P0)**2)/sum(P0SG*size) + sum(size*(P1SG-P1)**2)/sum(P1SG*size))
      pow[b,] <- sapply(alpha, function(z) pchisq(q=qchisq(z, df=length(size)-1, lower.tail=FALSE), ncp=ncp, df=length(size)-1, lower.tail=FALSE))
  }
  pow <- colMeans(pow) ; names(pow) <- alpha 
  pow
}  
  
