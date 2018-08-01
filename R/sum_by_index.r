sum_by_index <- function(X, INDEX) {
  if(!is.factor(INDEX)) stop("INDEX is not a factor")
  r <- .Call("oz_sum_by_group", PACKAGE = "Ravages", X, INDEX)
  names(r) <- levels(INDEX)
  r
}
