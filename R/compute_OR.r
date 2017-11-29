# Différents OR mais mêmes variants
OR.matrix <- function(n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro) {
  if(length(OR.del) != length(OR.pro))
    stop("Dimensions mismatch")
  OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
  # neutral, deleterious or protective
  v <- sample(1:3, n.variants, TRUE, c(1-prob.del-prob.pro, prob.del, prob.pro))
  t(apply(OR, 1, function(or) or[v]))
}
#example
#0R.matrix(20 , c(2,4), c(0.5,0.25), 0.2, 0.1)



# Mêmes OR mais variants différents
same.OR.matrix <- function(n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro){
  if(length(OR.del) != length(OR.pro))
    stop("Dimensions mismatch")
  OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
  v <- replicate(length(OR.del), sample(1:3, n.variants, TRUE, c(1 - prob.del - prob.pro, prob.del, prob.pro)))
  y <- matrix(rep(NA, n.variants*length(OR.del)), nrow=length(OR.del))
  for (i in 1:length(OR.del)) {y[i,] <- OR[i,][v[,i]]}
  return(y)
  }
#example
#same.OR.matrix(20, c(2,2), c(0.5,0.5), 0.2, 0.1)
