W.mat <- function(pi) {
  nind <- nrow(pi)
  ngpe <- ncol(pi)

  M <- vector("list", ngpe)
  for(i in 1:ngpe) M[[i]] <- do.call(rbind, lapply(1:ngpe, function(j) diag(-pi[,i]*pi[,j])))
  W <- do.call(cbind, M)

  diag(W) <- as.vector(pi*(1-pi))

  W
}


