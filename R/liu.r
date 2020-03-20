moments <- function(A, P){
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
  list(mu = muQ, sigma = sigmaQ, skewness = beta1, kurtosis = beta2)
}


p.valeur <- function(Q, moments, useskew) {
  s1 <- moments$skewness / sqrt(8)
  s2 <- moments$kurtosis / 12
  
  if(s1^2 > s2){
    a <- 1/(s1-sqrt(s1**2-s2))
    d <- s1*a**3-a**2
    l=a**2-2*d
  }else{
    if(useskew){
      a=1/s1
      d=0
      l=1/s1**2
      if(d <0 | l <0){ #Solution de secours
        l = 1/s2
        a = sqrt(l)
        d = 0
      }
    }else{
      l = 1/s2
      a = sqrt(l)
      d = 0
    }
  }
  muX <- l+d
  sigmaX <- sqrt(2)*a
  Q.norm <- (Q-moments$mu)/moments$sigma
  Q.norm1 <- Q.norm*sigmaX + muX
  
  #Calcul de la p-valeur suivant si l et d sont positifs
  #if(d>=0 & l>=0)  
  p.value <- pchisq(Q.norm1, df=l, ncp=d, lower.tail=FALSE)
  #else p.value <- p.valeur.withoutskew(Q = Q, moments = moments) ##Solution de secours
  return(p.value)
}

###Methode basée que sur les 3 moments
p.valeur.withoutskew.SKAT <- function(Q, moments){
  #Degre de liberte
  l <- ifelse(moments$kurtosis > 0, moments$kurtosis, 1e5)
  #Valeur du Chi-Deux
  chi2val <- l + sqrt(2 * l/moments$sigma) * (Q - moments$mu)
  #P-valeur associee
  p.value <- pchisq(chi2val, df=l, lower.tail=FALSE)
  return(p.value) 
}
