rbm.haplos.power <- function(haplos, freqs, weights = "SKAT", max.maf.causal = 0.01, maf.filter = max.maf.causal, p.causal = 0.5, p.protect = 0, h2 = c(0.01, 0.01), prev = c(1, 0.01), normal.approx = TRUE, size = c(500, 500), verbose = TRUE, alpha = 2.5e-6, RVAT = c("CAST", "WSS", "SKAT"), SKAT.method = c("permutations", "theoretical"), simus.haplos = c("freqs", "liability"), replicates = 1000, rep.by.causal = 50, cores = 10) {
  SKAT.method <- match.arg(SKAT.method)
  simus.haplos <- match.arg(simus.haplos)
  
  if(!("CAST" %in% RVAT | "WSS" %in% RVAT | "SKAT" %in% RVAT)) stop("Power calculations are only available for 'CAST', 'WSS' and 'SKAT'")
  #Simulation of data
  if(simus.haplos=="freqs"){
    x <- rbm.haplos.freqs(haplos = haplos, freqs = freqs, size = size, replicates = replicates)
  }
  if(simus.haplos=="liability"){
    x <- rbm.haplos.thresholds(haplos = haplos, weights = weights, max.maf.causal = max.maf.causal, p.causal = p.causal, p.protect = p.protect, h2 = h2, prev = prev, normal.approx = normal.approx, size = size, replicates = replicates, rep.by.causal =rep.by.causal, verbose = verbose)
  }
  #Filtering to keep only rare variants
  x <- filter.rare.variants(x, maf.threshold = maf.filter, min.nb.snps = 2)
  ###RVAT
  pow.names <- c()
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
      x.WSS.pow <- sapply(alpha, function(z) mean(x.WSS$p.value<z, na.rm = T) )
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
    x.SKAT.pow <- sapply(alpha, function(z) mean(x.SKAT$p.value<z))
    pow.names <- c(pow.names, "SKAT")
  }else{
    x.SKAT.pow <- NULL
  }
  pow <- rbind(x.CAST.pow, x.WSS.pow, x.SKAT.pow)
  rownames(pow) <- pow.names
  colnames(pow) <- alpha
  pow
}
