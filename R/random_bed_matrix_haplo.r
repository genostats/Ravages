
random.bed.matrix.haplos <- function(haplos, freqs, size, B) {
  bed <- .Call('random_bed_matrix_haplotypes_freqs', PACKAGE = "Ravages", haplos, freqs, size, B)

  nb_inds <- sum(size);

  ids <- sprintf("A%0*d", log10(nb_inds) + 1, 1:nb_inds)
  ped <- data.frame(famid = ids,  id = ids, father = 0, mother = 0, sex = 0,
            pheno = unlist(mapply(rep, 1:length(size), each = size, SIMPLIFY = FALSE)) - 1, 
            stringsAsFactors = FALSE)

  snps <- data.frame(chr = NA, id = NA, dist = NA, pos = NA, A1 = NA, A2 = NA, 
               genomic.region = factor( rep(sprintf("R%0*d", log10(B) + 1, 1:B), each = ncol(haplos)) ),
               stringsAsFactors = FALSE)
 
  if( is.null(colnames(haplos)) )
    snps$id <- paste( snps$genomic.region, sprintf("SNP%0*d", log10(ncol(haplos)) + 1, 1:ncol(haplos)), sep = "_")
  else
    snps$id <- paste( snps$genomic.region, colnames(haplos), sep = "_")

  x <- new("bed.matrix", bed = bed, snps = snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x, verbose = FALSE)
  x
}

random.bed.matrix.haplos.thresholds <- function(haplos, burdens, sd, thr1, thr2, size, B) {
  bed <- .Call('random_bed_matrix_haplotypes_thresholds', PACKAGE = "Ravages", haplos, burdens, sd, thr1, thr2, size, B)

  nb_inds <- sum(size);

  ids <- sprintf("A%0*d", log10(nb_inds) + 1, 1:nb_inds)
  ped <- data.frame(famid = ids,  id = ids, father = 0, mother = 0, sex = 0,
            pheno = unlist(mapply(rep, 1:length(size), each = size, SIMPLIFY = FALSE)) - 1, 
            stringsAsFactors = FALSE)

  snps <- data.frame(chr = NA, id = NA, dist = NA, pos = NA, A1 = NA, A2 = NA, 
               genomic.region = factor( rep(sprintf("R%0*d", log10(B) + 1, 1:B), each = ncol(haplos)) ),
               stringsAsFactors = FALSE)
 
  if( is.null(colnames(haplos)) )
    snps$id <- paste( snps$genomic.region, sprintf("SNP%0*d", log10(ncol(haplos)) + 1, 1:ncol(haplos)), sep = "_")
  else
    snps$id <- paste( snps$genomic.region, colnames(haplos), sep = "_")

  x <- new("bed.matrix", bed = bed, snps = snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x, verbose = FALSE)
  x
}

