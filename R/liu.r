p.valeur <- function(A, P, Q){
  A1 <- A %*% P
  # on pourrait utiliser t(G) %*% P %*% G  [ si A = GG' ]
  A2 <- A1 %*% A1
  # les traces [on a delta_i = 0]
  c1 <- c(sum(diag(A1)), sum(diag(A2)), sum(diag(A1 %*% A2)), sum(diag(A2 %*% A2)))
  
  muQ <- c1[1]
  sigmaQ <- sqrt(2*c1[2])
  s1 <- c1[3]/c1[2]^(3/2)
  s2 <- c1[4]/c1[2]^2
  beta1 <- sqrt(8) * s1
  beta2 <- 12*s2
 
  if(s1^2 > s2){
    a <- 1/(s1-sqrt(s1**2-s2))
    d <- s1*a**3-a**2
    l=a**2-2*d
  }else{
    a=1/s1
    d=0
    l=1/s1**2
  }
  muX <- l+d
  sigmaX <- sqrt(2)*a
  Q.norm <- (Q-muQ)/sigmaQ
  Q.norm1 <- Q.norm*sigmaX + muX
  
  #Calcul de la p-valeur
  p.value <- pchisq(Q.norm1, df=l, ncp=d, lower.tail=FALSE)
  return(p.value)
}
