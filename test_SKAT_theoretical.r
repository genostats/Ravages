require(Ravages)
example(rbm.haplos.freqs)

if(FALSE) {
set.seed(1)
x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size = rep(10, 5), replicates = 10)
NO <- NullObject.parameters(x@ped$pheno, RVAT = "SKAT", pheno.type = "categorial")
system.time(a <- SKAT.theoretical(x, NO, cores = 1))
system.time(b <- Ravages:::SKAT.theoretical1(x, NO, cores = 1))
}


set.seed(1)
x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size = rep(100, 5), replicates = 10) 
NO <- NullObject.parameters(x@ped$pheno, RVAT = "SKAT", pheno.type = "categorial")
system.time(a <- SKAT.theoretical(x, NO, cores = 1, debug = TRUE))
system.time(b <- Ravages:::SKAT.theoretical1(x, NO, cores = 1, debug = TRUE))

# 26.199   0.664  26.883 
#  7.069   0.112   7.159 
 
if(FALSE){
set.seed(1)
x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size = rep(200, 5), replicates = 10) 
NO <- NullObject.parameters(x@ped$pheno, RVAT = "SKAT", pheno.type = "categorial")
system.time(a <- SKAT.theoretical(x, NO, cores = 1))
system.time(b <- Ravages:::SKAT.theoretical1(x, NO, cores = 1))
}
