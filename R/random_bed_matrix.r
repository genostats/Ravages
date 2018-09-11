
# il faut ajouter un argument pour la fonction à appeler pour générer les OR
# (ici OR.matrix...)
random.bed.matrix <- function(pop.maf, size, baseline, replicates, OR.pars, OR.function) {
  OR.pars$n.variants <- length(pop.maf)
  nb_snps <- OR.pars$n.variants * replicates
  nb_inds <- sum(size)
  x <- new.bed.matrix(nb_inds, nb_snps);
  for(b in 1:replicates) {
    OR <- do.call( OR.function, OR.pars)
    MAFS <- group.mafs(pop.maf, OR, baseline)
    .Call("oz_random_filling_bed_matrix", PACKAGE = "Ravages", x@bed, MAFS, size, (b-1)*OR.pars$n.variants)
  }
  x@ped$pheno <- rep.int( 1:length(size) - 1, size)
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = OR.pars$n.variants) )
  x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  x <- set.stats(x, verbose = FALSE)
  x
}

