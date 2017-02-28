
burden <- function(x, regions, weights.0, weights.1, weights.2) {
  if(length(weights.1) != ncol(x) | length(weights.2) != ncol(x) | length(regions) != ncol(x)) {
    stop("x and weights dimensions mismatch")
  }
  if(!is.factor(regions)) stop("'regions' is not a factor")

  B <- .Call('oz_burden2', PACKAGE = 'oz', x@bed, nlevels(regions), regions, weights.0, weights.1, weights.2)
  colnames(B) <- levels(regions)
  rownames(B) <- x@ped$id
  B
}
