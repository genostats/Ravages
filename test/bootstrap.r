source("~/working_directory/Ravages/test/grandes_matrices.r")

# ymp = y - pi^

Qstat <- function(ymp, K, n) {
  a <- apply(ymp, 2, function(x) x %*% K %*% x) 
  a[n == 0] <- NA
  sum(a/n, na.rm = TRUE)
}

# "parametric bootstrap" (sort of ?)
# -> generation de valeurs de Q à parti de la matrice pi (telle que la sort Pi.matrix)
# si on fournit n (n = table(Y)) on gardera ces valeurs fixées, sinon c'est recalculé à chaque tirage.
bootstrap <- function(B, pi, X, K, n) {
  recompute_n <- if(missing(n)) TRUE else FALSE
  W <- W.mat(pi[,-1,drop=FALSE])
  XX <- block.diag( rep(list(X), ncol(pi) - 1) )
  WX <- W %*% XX
  XWX <- t(XX) %*% WX
  U <- WX %*% solve(XWX, t(XX))

  lev <- colnames(pi)
  if(is.null(lev))
    lev = sprintf("L%d", seq(0, length = ncol(pi)))

  R <- numeric(B)
  for(b in 1:B) {
    Y <- factor( apply(pi, 1, function(p) sample(lev, 1, prob = p)), levels = lev)

    if(recompute_n)
      n <- table(Y)
   
    # Y - pi^_0 où pi^_0  est le vecteur des pi^ estimé sur les données
    YY <- sapply(lev, function(l) as.numeric(Y == l))
    ymp0 <- as.vector(YY[,-1]) - as.vector(pi[,-1])
    # et Y - pi^ "bootstrapé"
    ymp <- matrix( ymp0 - U %*% ymp0, nrow = nrow(pi))
    ymp <- cbind( -rowSums(ymp), ymp ) # la composante pour le groupe 0
    # finalement, la stat
    R[b] <- Qstat(ymp, K, n)
  }
  
  m <- mean(R)
  s <- sd(R)
  Rs <- (R - m)/s
  beta1 <- mean( Rs**3 )
  beta2 <- mean( Rs**4 ) - 3
  return(list( bootstrap = R, mean = m, sigma = s, skewness = beta1, kurtosis = beta2 ) )
}


# test des fonctions cpp
Rcpp::sourceCpp("~/working_directory/Ravages/test/bootstrap.cpp")
bootstrap_cpp <- function(B, pi, X, K) {
  W <- W.mat(pi[,-1,drop=FALSE])
  XX <- block.diag( rep(list(X), ncol(pi) - 1) )
  WX <- W %*% XX
  XWX <- t(XX) %*% WX
  U <- -(WX %*% solve(XWX, t(XX)))
  diag(U) <- diag(U)+1

  lev <- colnames(pi)
  if(is.null(lev))
    lev = sprintf("L%d", seq(0, length = ncol(pi)))

  R <- numeric(B)
  for(b in 1:B) {
    r <- bootstrap1(U, pi)
    R[b] <- Qstat(r$ymp, K, r$n)
  }
  
  m <- mean(R)
  s <- sd(R)
  Rs <- (R - m)/s
  beta1 <- mean( Rs**3 )
  beta2 <- mean( Rs**4 ) - 3
  return(list( bootstrap = R, mean = m, sigma = s, skewness = beta1, kurtosis = beta2 ) )
}
