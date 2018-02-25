require(oz)
my.pars <- list( OR.del = c(2,4), prob.del = 0.2, prob.pro = 0.1)
oz:::fff( c(0.01, 0.02, 0.04), c(200, 100, 100), c(0.01, 0.01), 10, my.pars ) -> x

CAST(x, x@snps$genomic.region, group = x@ped$pheno) 
WSS(x)
C.ALPHA(x, B.max=100)
Beta.M(x, B.max=100)

Power( pop.maf = c(0.01, 0.02, 0.04), size = c(200, 100, 100), baseline = c(0.01, 0.01), replicates = 100, OR.pars = my.pars ) 


# ------
require(oz)
my.OR.pars <- list(OR.del=c(3,6), prob.del=0.2, prob.pro=0.05)
x <- random.bed.matrix(Kryukov$maf[1:50], c(400,200,200), c(0.001,0.001), 10, OR.pars=my.OR.pars, scenario=2)
Beta.M(x)
Beta.M.rect(x)

 p1 <- Beta.M(x, target=500)$p.value
 p2 <- Beta.M.rect(x, target=500)$p.value


 x <- random.bed.matrix(Kryukov$maf[1:50], c(400,200,200), c(0.001,0.001), 100, OR.pars=my.OR.pars, scenario=1)
 s1 <- Beta.M(x)$stat
 s2 <- Beta.M.rect(x)$stat
 plot(s1,s2)

# ------
require(oz)
my.OR.pars <- list(OR.del=c(3,6), prob.del=0.2, prob.pro=0.05)
x <- random.bed.matrix(Kryukov$maf[1:50], c(400,200,200), c(0.001,0.001), 10, OR.pars=my.OR.pars, scenario=2)
Beta.M(x)
Sum.Fst(x)

 p1 <- Beta.M(x, target=500)$p.value
 p2 <- Sum.Fst(x, target=500)$p.value 

