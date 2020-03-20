block.diag <- function(L) {
  if(length(L) == 1)
    return(L[[1]])
  D <- block.diag( L[-1] )
  L1 <- L[[1]]
  ZerosNE <- matrix(0, nrow = nrow(L1), ncol = ncol(D))
  ZerosSW <- matrix(0, nrow = nrow(D), ncol(L1))
  rbind( cbind(L1, ZerosNE),   cbind(ZerosSW, D) )
}
