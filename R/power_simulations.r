
# un des scenario d'OR ...
OR.matrix <- function(n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro) {
  if(length(prob.del) != length(prob.pro)) 
    stop("Dimensions mismatch")
  OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
  # neutral, deleterious or protective
  v <- sample(1:3, n.variants, TRUE, c(1-prob.del-prob.pro, prob.del, prob.pro))
  t(apply(OR, 1, function(or) or[v]))
}
#example
#0R.matrix(20 , c(2,4), c(0.5,0.25), 0.2, 0.1)

# pas inspiré pour trouver un nom
# et y a trop de paramètres c'est un peu le bordel
fff <- function(pop.maf, size, baseline, replicates, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro) {
  ff <- function() {
    OR <- OR.matrix( length(pop.maf), OR.del, OR.pro, prob.del, prob.pro )
    MAFS <- group.mafs(pop.maf, OR, baseline)
    random.bed.matrix(MAFS, size)
  }
  x <- suppressWarnings(do.call( cbind, replicate(replicates, ff())))  # warnings à cause 
  x@snps$genomic.region <- rep(sprintf("R%d", 1:replicates), each = length(pop.maf))
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x
}


