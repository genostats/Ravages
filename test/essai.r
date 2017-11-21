require(oz)
my.pars <- list( OR.del = c(2,4), prob.del = 0.2, prob.pro = 0.1)
oz:::fff( c(0.01, 0.02, 0.04), c(200, 100, 100), c(0.01, 0.01), 10, my.pars ) -> x

CAST(x, x@snps$genomic.region, group = x@ped$pheno) 
WSS(x)
C.ALPHA(x, B.max=100)
Beta.M(x, B.max=100)

Power( pop.maf = c(0.01, 0.02, 0.04), size = c(200, 100, 100), baseline = c(0.01, 0.01), replicates = 100, OR.pars = my.pars ) 

