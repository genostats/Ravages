
random.bed.matrix.genomic.probs <- function(pop.proba0, pop.proba1, size, replicates) {
  nb_variants <- ncol(pop.proba0) 
  nb_snps <- nb_variants * replicates
  nb_inds <- sum(size)
  x <- new.bed.matrix(nb_inds, nb_snps);
  for(b in 1:replicates) {
    .Call("oz_random_filling_bed_matrix_noHW", PACKAGE = "Ravages", x@bed, pop.proba0, pop.proba1, size, (b-1)*nb_variants)
  }
  x@ped$pheno <- rep.int( 1:length(size) - 1, size)
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = nb_variants) )
  x@snps$id <- paste( x@snps$genomic.region, sprintf("m%d", 1:nb_variants), sep="_")
  x <- set.stats(x, verbose = FALSE)
  x
}

