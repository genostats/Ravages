require(Ravages)
example(rbm.haplos.freqs)

x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size = rep(10, 5), replicates = 10)
NO <- NullObject.parameters(x@ped$pheno, RVAT = "SKAT", pheno.type = "categorial")
system.time(a <- SKAT.theoretical(x, NO, cores = 1))
system.time(b <- Ravages:::SKAT.theoretical1(x, NO, cores = 1))

stop()
x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size = rep(100, 5), replicates = 10) 
NO <- NullObject.parameters(x@ped$pheno, RVAT = "SKAT", pheno.type = "categorial")
system.time(a <- SKAT.theoretical(x, NO, cores = 1))
system.time(b <- Ravages:::SKAT.theoretical1(x, NO, cores = 1))

x <- rbm.haplos.freqs(haplos=LCT.gene.hap, freqs=LCT.freqs, size = rep(200, 5), replicates = 10) 
NO <- NullObject.parameters(x@ped$pheno, RVAT = "SKAT", pheno.type = "categorial")
system.time(a <- SKAT.theoretical(x, NO, cores = 1))
system.time(b <- Ravages:::SKAT.theoretical1(x, NO, cores = 1))
