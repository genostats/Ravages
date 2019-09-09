simulated.bedmatrix.haplo <- function(haplos, weights.variants=weights, maf.threshold = 0.01, p.causal, p.protect = 0, nb.causal, h2, prev, normal.approx = TRUE, scenario, size = c(1000, 500, 500), replicates = 10, power, alpha.threshold=2.5e-6){
  if(!is.matrix(haplos)) stop("haplos is not a matrix")

  if(!is.function(weights.variants)){
    stop("weighths.variants is not a function")
  }else{
    weights.var <- do.call(weights.variants, list(haplos, maf.threshold))
  }

  if(!(scenario %in% c("SVDR", "DG", "SVSR", "DVSG"))) stop(paste("scenario", scenario, "doesn't exist"))

  #Get burdens, thresholds, ...
  if(scenario != "DVSG"){
    RR <- simu.prolego(haplos, weights.var, nb.causal = nb.causal, h2 = h2, prev = prev) #1 seul groupe d'individus
    x <- do.call(scenario, list(RR, size, replicates))

  }else{
    RR1 <- simu.prolego(haplos, weights.var, nb.causal = nb.causal, h2 = h2, prev = prev)
    RR2 <- simu.prolego(haplos, weights.var, nb.causal = nb.causal, h2 = h2, prev = prev)
    x <- do.call(scenario, list(RR1, RR2, size, replicates))
  }

  if(missing(power)) return(x)

  if(!missing(power)){
    if(!(power %in% c("CAST", "WSS", "SKAT"))) stop("Power calculation only for CAST, WSS or SKAT")

    #Power only on rare variants
    x <- filter.rare.variants(x, maf.threshold=maf.threshold, filter="whole")

    power.res <- data.frame(ThreeG=NA, TwoG=NA, Group1=NA, Group2=NA, test=rep(power, each=replicates))
    if("CAST" %in% power){ power.res[power.res$test=="CAST", 1:4] <- run.CAST(x) }
    if("WSS" %in% power){ power.res[power.res$test=="WSS", 1:4] <- run.WSS(x) }
    if("SKAT" %in% power){ power.res[power.res$test=="SKAT", 1:4] <- run.skat(x) }
    
    pow.values <- sapply(unique(power.res$test), function(tt) apply(subset(power.res, test==tt)[,1:4], 2, function(analysis) mean(analysis<alpha.threshold, na.rm=TRUE)))
    colnames(pow.values) <- unique(power)
    return(list(power=pow.values,x=x))
  }

}
