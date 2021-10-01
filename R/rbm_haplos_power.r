rbm.haplos.power <- function(haplos, freqs, weights = -0.4*log10(colMeans(haplos)), maf.threshold = 0.01, nb.causal, p.protect = 0, h2, prev, normal.approx = TRUE, size, verbose = TRUE, alpha = 2.5e-6, RVAT = c("CAST", "WSS", "SKAT"), simus.haplos = c("freqs", "liability"), replicates = 1000, rep.by.causal = 50, cores = 10) {
  simus.haplos <- match.arg(simus.haplos)
  if(!("CAST" %in% RVAT | "WSS" %in% RVAT | "SKAT" %in% RVAT)) stop("Power calculations are only available for 'CAST', 'WSS' and 'SKAT'")
  #Simulation of data
  if(simus.haplos=="freqs"){
    x <- rbm.haplos.freqs(haplos = haplos, freqs = freqs, size = size, replicates = replicates)
  }
  if(simus.haplos=="liability"){
    x <- rbm.haplos.thresholds(haplos = haplos, weights = weights, maf.threshold = maf.threshold, nb.causal = nb.causal, p.protect = p.protect, h2 = h2, prev = prev, normal.approx = normal.approx, size = size, replicates = replicates, rep.by.causal = 50, verbose = verbose)
  }
  #Keeping only rare variants
  x <- filter.rare.variants(x, filter = "whole", maf.threshold = maf.threshold)
  ###RVAT
  pow.names <- c()
  if("CAST" %in% RVAT | "WSS" %in% RVAT){
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
  pow <- c(x.CAST.pow, x.WSS.pow, x.SKAT.pow)
  names(pow) <- pow.names
  pow
}