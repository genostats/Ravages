require(Ravages)


     x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)
     #Add population
     x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]
     
     #Select EUR superpopulation
     x <- select.inds(x, superpop=="EUR")
  
     w <- c( head(which(x@ped$pop == "CEU"), 10), head(which(x@ped$pop == "TSI"), 5), head(which(x@ped$pop == "FIN"), 20) )
     x <- x[w,]
     x@ped$pop <- droplevels(x@ped$pop)

     
     #Group variants within known genes
     x <- set.genomic.region(x)
     
     #Filter of rare variants: only non-monomorphic variants with
     #a MAF lower than 2.5%
     #keeping only genomic regions with at least 10 SNPs
     x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)

if(FALSE) {
     sex <- x1@ped$sex
     set.seed(1) ; u <- runif(nrow(x1))
     D <- cbind(sex, u)

     x1.H0 <- NullObject.parameters(x1@ped$pop, RVAT = "SKAT", pheno.type = "categorial", data = D)

     debug(Ravages:::get.parameters.pvalue.theoretical)
     SKAT(x1, x1.H0, get.moments = "theoretical", cores = 1)
}

# à partir de ça on a récupéré P = x1.H0$P1, we = les poids des SNPs, et A = G'PG avec G = (la région R3HDM1 avec poids we)
# Browse[2]> GPG <<- ( t(G.bloc) %*% P1 %*% G.bloc )
# Browse[2]> we <<- x.genomic.region@snps$weights # (avant l'extraction)
# Browse[2]> P <<- P1
# save(GPG, we, P, file = "test_blocs.rda")
load("test_blocs.rda")
Z <- .Call("Ravages_GPG", PACKAGE = "Ravages", x1@bed, (x1@snps$genomic.region == "R3HDM1"), x1@p, we, P, c(10L, 5L, 20L))
range(Z - GPG)

