random.bed.matrix.haplos <- function(haplos, weights = burden.weights, maf.threshold = 0.01, p.causal, p.protect = 0, nb.causal, h2, prev, 
                                      size, replicates = 1){
  if(!is.matrix(haplos)) stop("haplos is not a matrix")

  if(!is.function(weights)){
    stop("weights is not a function")
  }else{
    weights.var <- do.call(weights, list(haplos, maf.threshold))
  }

  if(!missing(p.causal) & !missing(nb.causal))
    stop("Please give either nb.causal or p.causal but not both")

  # on genere un vecteur de directions d'effet a partir de p.causal et p.protect
  # (sampling qui impose exactement les proportions demandees) !
  w <- which(weights.var > 0) # parmi les SNPs potentiellement causaux
  if(missing(nb.causal)){
    nb.causal <- ceiling(p.causal*length(w))
  }

  n.group <- max( length(nb.causal), length(p.protect), length(h2), length(prev), length(size) )
  if( any( !(c( length(nb.causal), length(p.protect), length(h2), length(prev), length(size) ) %in% c(1, n.group) ) ) ) 
    stop("Arguments length mismatch")
  if(length(nb.causal) < n.group) nb.causal <- rep(nb.causal, n.group)
  if(length(p.protect) < n.group) p.protect <- rep(p.protect, n.group)
  if(length(h2) < n.group) h2 <- rep(h2, n.group)
  if(length(prev) < n.group) prev <- rep(prev, n.group)
  if(length(size) < n.group) size <- rep(size, n.group)

  RR <- mapply(simu.prolego, list(haplos), list(weights.var), p.protect, nb.causal, h2, prev, TRUE, SIMPLIFY = FALSE)
  RR <- list(haplos = haplos, burdens = sapply(RR, function(x) x$burdens), thr1 = sapply(RR, function(x) x$thr1), thr2 = sapply(RR, function(x) x$thr2), 
             sd = sqrt(1-h2), size = size, B = replicates)
  x <- do.call(rbm.haplos.thresholds, RR)
}
