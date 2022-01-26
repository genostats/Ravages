
GPG <- function(x, we, P, N = c(10L, 5L, 20L)) {
  .Call("Ravages_GPG", PACKAGE = "Ravages", x@bed, (x@snps$genomic.region == "R3HDM1"), x@p, we, P, N) 
}

