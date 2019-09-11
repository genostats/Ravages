simulated.bedmatrix.haplo <- function(haplos, weights.variants=burden.weights, maf.threshold = 0.01, p.causal, p.protect = 0, nb.causal, h2, prev, normal.approx = TRUE, scenario, size = c(1000, 500, 500), replicates = 10, power, alpha.threshold=2.5e-6){
  if(!is.matrix(haplos)) stop("haplos is not a matrix")

  if(!is.function(weights.variants)){
    stop("weighths.variants is not a function")
  }else{
    weights.var <- do.call(weights.variants, list(haplos, maf.threshold))
  }

  if(!(scenario %in% c("SVDR", "DG", "SVSR", "DVDR", "DVSR"))) stop(paste("scenario", scenario, "doesn't exist"))
  
  #Get burdens, thresholds, ...
  if(scenario != "DVDR"){
    if(length(prev)>1 | length(h2)>1){
      if(scenario == "DG"){
        warning("Only second group of cases with genetic effects: only one prevalence and h2 needed")
      }else{
        warning("Same prevalence and h2 value between the groups of cases in this scenario, only the first one is used")
      }
      prev <- prev[1] ; h2 <- h2[1]
    }
    RR <- simu.prolego(haplos, weights.var, nb.causal = nb.causal, h2 = h2, prev = prev) #1 seul groupe d'individus
    x <- do.call(scenario, list(RR, size, replicates))

  }else{
    if(length(prev)==1 | length(h2)==1 | length(prev) != length(h2)) stop("Needs one prevalence and one h2 per group of cases")
    #RR1 <- simu.prolego(haplos, weights.var, nb.causal = nb.causal, h2 = h2[1], prev = prev[1])
    #RR2 <- simu.prolego(haplos, weights.var, nb.causal = nb.causal, h2 = h2[2], prev = prev[2])
    #x <- do.call(scenario, list(RR1, RR2, size, replicates))
    RR <- lapply(1:length(prev), function(z) simu.prolego(haplos, weights.var, nb.causal = nb.causal, h2 = h2[z], prev=prev[z]))
    names(RR) <- paste("RR", 1:length(prev), sep="")
    x <- do.call(scenario, list(RR, size, replicates))
  }

  if(missing(power)) return(x)

  if(!missing(power)){
    if(length(intersect(power, c("CAST", "WSS", "SKAT")))<length(power)) stop("Power calculation only for CAST, WSS or SKAT")

    #Power only on rare variants
    x <- filter.rare.variants(x, maf.threshold=maf.threshold, filter="whole")

    power.res <- data.frame(ThreeG=NA, TwoG=NA, Group1=NA, Group2=NA, test=rep(power, each=replicates))
    if("CAST" %in% power){ power.res[power.res$test=="CAST", 1:4] <- run.CAST(x) }
    if("WSS" %in% power){ power.res[power.res$test=="WSS", 1:4] <- run.WSS(x) }
    if("SKAT" %in% power){ power.res[power.res$test=="SKAT", 1:4] <- run.skat(x) }
    
    pow.values <- sapply(unique(power.res$test), function(tt) apply(subset(power.res, power.res$test==tt)[,1:4], 2, function(analysis) mean(analysis<alpha.threshold, na.rm=TRUE)))
    colnames(pow.values) <- unique(power)
    return(list(power=pow.values,x=x))
  }

}
