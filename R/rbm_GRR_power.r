rbm.GRR.power <- function(genes.maf = Kryukov, size, prev, GRR.matrix.del, GRR.matrix.pro, p.causal = 0.5, p.protect = 0, same.variant = FALSE, genetic.model = c("general", "multiplicative", "dominant", "recessive"), select.gene, alpha = 2.5e-6, selected.controls = TRUE, power.type = c("simulations", "theoretical"), verbose = TRUE, RVAT = c("CAST", "WSS", "SKAT"), maf.threshold = 0.01, replicates = 1000, cores = 10){
  power.type <- match.arg(power.type)
  if(power.type != "simulations" & power.type != "theoretical") stop("'power.type' should be 'simulations' or 'theoretical'")
  if(power.type=="simulations" & !("CAST" %in% RVAT | "WSS" %in% RVAT | "SKAT" %in% RVAT)) stop("Power calculations are only available for 'CAST', 'WSS' and 'SKAT'")
  if (nlevels(genes.maf$gene) > 1) {
    if (missing(select.gene)) {
      warning("More than one gene in the file, only the first one is used")
      select.gene <- levels(genes.maf$gene)[[1]]
    }
    pop.maf <- subset(genes.maf, genes.maf$gene %in% select.gene)$maf
  } else {
    pop.maf <- genes.maf$maf
  }
  if (!is.list(GRR.matrix.del)) {
    if (is.matrix(GRR.matrix.del)) {
      GRR.matrix.del <- list(GRR.matrix.del)
    } else {
      stop("GRR.matrix.del should be a list or a matrix")
    }
    GRR.het <- GRR.matrix.del[[1]]
  }
  genetic.model <- match.arg(genetic.model)
  if (length(GRR.matrix.del) == 1) {
    if (genetic.model == "general") {
      stop("Needs two GRR matrices in the general model")
    } else {
      GRR.homo.alt <- NULL
    }
  } else {
    if (genetic.model == "general") {
      GRR.homo.alt <- GRR.matrix.del[[2]]
    } else {
      warning("Only one GRR matrix needed for this model, only the first one is used")
      GRR.homo.alt <- NULL
    }
  }
  if (!missing(GRR.matrix.pro)) {
    if (!is.list(GRR.matrix.pro)) {
      if (is.matrix(GRR.matrix.pro)) {
        GRR.matrix.pro <- list(GRR.matrix.pro)
      } else {
        stop("GRR.matrix.pro should be a list or a matrix")
      }
    }
    GRR.het.pro <- GRR.matrix.pro[[1]]
    if (length(GRR.matrix.pro) == 1) {
      if (genetic.model == "general") {
        stop("Needs two GRR matrices in the general model")
      } else {
        GRR.homo.alt.pro <- NULL
      }
    } else {
      if (genetic.model == "general") {
        GRR.homo.alt.pro <- GRR.matrix.pro[[2]]
        } else {
          warning("Only one GRR matrix needed for this model, only the first one is used")
          GRR.homo.alt.pro <- NULL
        }
      }
  } else {
    GRR.het.pro <- 1/GRR.het
    GRR.homo.alt.pro <- 1/GRR.homo.alt
  }
  if (any(GRR.het < 1) | any(GRR.homo.alt < 1)) stop("Matrix of deleterious GRR has GRR values lower than 1")
  if (!missing(GRR.matrix.pro)) {
    if (any(GRR.het.pro > 1) | any(GRR.homo.alt.pro > 1))
      stop("Matrix of protective GRR has GRR values greater than 1")
  }
  GRR.pars <- list(OR.del = GRR.het, OR.pro = GRR.het.pro, p.causal = p.causal, prob.pro = p.protect)
  if (!is.null(GRR.homo.alt)) {
    GRR.homo.alt.pars <- list(OR.del = GRR.homo.alt, OR.pro = GRR.homo.alt.pro)
  } else {
    GRR.homo.alt.pars <- NULL
  }
  if (!is.null(GRR.homo.alt.pars)) {
    if (nrow(GRR.het) != nrow(GRR.homo.alt.pars$OR.del) |
      ncol(GRR.het) != ncol(GRR.homo.alt.pars$OR.del))
      stop("GRR.het and GRR.homo.alt have different dimensions")
  }
  if (same.variant) {
    variant.function <- OR.matrix.same.variant
  } else {
    variant.function <- OR.matrix
  }
  GRR.pars$n.variants <- length(pop.maf)
  nb_snps <- GRR.pars$n.variants * replicates
  nb_inds <- sum(size)

  if(power.type=="simulations") x <- new.bed.matrix(nb_inds, nb_snps)
  else pow <- c()
  for (b in 1:replicates) {
      GRR.het <- do.call(variant.function, GRR.pars)
      if (!is.null(GRR.homo.alt.pars)) {
        for (i in 1:length(prev)) {
          GRR.homo.alt[i, which(GRR.het[i, ] == 1)] <- 1
          GRR.homo.alt[i, which(GRR.het[i, ] < 1)] <- GRR.homo.alt.pars$OR.pro[i, which(GRR.het[i, ] < 1)]
        }
      } else {
        if(power.type=="simulations") GRR.homo.alt = NULL
        else GRR.homo.alt = GRR.het**2
      }
    if(power.type=="theoretical"){
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
      ##Compute power 
      pow <- c(pow, pchisq(q=qchisq(alpha, df=length(size)-1, lower.tail=FALSE), ncp=ncp, df=length(size)-1, lower.tail=FALSE))
    }else{
    #Simulations as in rbm.GRR()
      MAFS <- genotypic.freq(genes.maf = genes.maf, GRR.het = GRR.het, GRR.homo.alt = GRR.homo.alt, prev = prev, select.gene = select.gene, genetic.model = genetic.model, selected.controls = selected.controls)
      if (any(MAFS$freq.homo.ref[1, ] > 1 | MAFS$freq.het[1,] < 0 | MAFS$freq.homo.alt[1, ] < 0)) stop("Impossible genetic model, please change your parametrization")
      .Call("oz_random_filling_bed_matrix_noHW", PACKAGE = "Ravages", x@bed, MAFS$freq.homo.ref, MAFS$freq.het, size, (b - 1) * GRR.pars$n.variants)
    }
  }
  if(power.type=="simulations"){
    x@ped$pheno <- factor(rep.int(1:length(size) - 1, size))
    x@snps$genomic.region <- factor(rep(sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = GRR.pars$n.variants))
    x@snps$id <- paste(x@snps$genomic.region, x@snps$id, sep = "_")
    #Keeping only rare variants
    x <- filter.rare.variants(x, filter = "whole", maf.threshold = maf.threshold)
    ##RVAT
    if("CAST" %in% RVAT | "WSS" %in% RVAT){
      pow.names <- c()
      H0.burden <- NullObject.parameters(x@ped$pheno, RVAT = "burden", pheno.type = "cat", ref.level = 0)
      if("CAST" %in% RVAT){
        x.CAST <- burden(x, H0.burden, burden = "CAST", verbose = verbose, maf.threshold = maf.threshold, cores = cores)
        x.CAST.pow <- mean(x.CAST$p.value<alpha, na.rm = T)
        pow.names <- c(pow.names, "CAST")
      }else{
        x.CAST.pow <- NULL
      }
      if("WSS" %in% RVAT){
        x.WSS <- burden(x, H0.burden, burden = "WSS", verbose = verbose, cores = cores)
        x.WSS.pow <- mean(x.WSS$p.value<alpha, na.rm = T) 
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
      x.SKAT <- SKAT(x, H0.SKAT, verbose = verbose, cores = cores)
      x.SKAT.pow <- mean(x.SKAT$p.value<alpha)
      pow.names <- c(pow.names, "SKAT")
    }else{
      x.SKAT.pow <- NULL
    }
  }
  if(power.type=="theoretical") pow <- mean(pow)
  else{ pow <- c(x.CAST.pow, x.WSS.pow, x.SKAT.pow) ; names(pow) <- pow.names }
  return(pow)
}
  
  
