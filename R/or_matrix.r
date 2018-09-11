# Mêmes OR mais variants différents
OR.matrix <- function(n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro){
  if(length(OR.del) != length(OR.pro))
    stop("Dimensions mismatch")
  OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
  v <- replicate(length(OR.del), sample(1:3, n.variants, TRUE, c(1 - prob.del - prob.pro, prob.del, prob.pro)))
  t(sapply(1:nrow(OR), function(x) OR[x,][v[,x]]))
  }
#example
#same.OR.matrix(20, c(2,2), c(0.5,0.5), 0.2, 0.1)
