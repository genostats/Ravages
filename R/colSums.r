

colSumsSq  <- function(x) .Call("oz_colsums_sq",  PACKAGE="oz", x)
colSumsCub <- function(x) .Call("oz_colsums_cub", PACKAGE="oz", x)

