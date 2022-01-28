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


     sex <- x1@ped$sex
     set.seed(1) ; u <- runif(nrow(x1))
     D <- cbind(sex, u)

     x1.H0 <- NullObject.parameters(x1@ped$pop, RVAT = "SKAT", pheno.type = "categorial", data = D)

SKAT.theoretical(x1, x1.H0, cores = 1, debug = TRUE)
Ravages:::SKAT.theoretical1(x1, x1.H0, cores = 1, debug = TRUE)

