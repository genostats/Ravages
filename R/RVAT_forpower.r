run.skat <- function(x) {

  skat <- Ravages:::SKAT(x, weights = SKAT:::Beta.Weights(x@snps$maf, c(1,25)), maf.threshold = 1)
  pool <- Ravages:::SKAT(x, ifelse(x@ped$pheno == 0, 0, 1), weights = SKAT:::Beta.Weights(x@snps$maf, c(1,25)), maf.threshold = 1)

  x.1 <- select.inds(x, x@ped$pheno != 2)
  x.2 <- select.inds(x, x@ped$pheno != 1)
  x.2@ped$pheno <- ifelse(x.2@ped$pheno == 0, 0, 1)
  cas1 <- Ravages:::SKAT(x.1, weights = SKAT:::Beta.Weights(x@snps$maf, c(1,25)), maf.threshold = 1)
  cas2 <- Ravages:::SKAT(x.2, weights = SKAT:::Beta.Weights(x@snps$maf, c(1,25)), maf.threshold = 1)

  data.frame(skat = skat$p.value, skat.pool = pool$p.value, skat.1 = cas1$p.value, skat.2 = cas2$p.value)
}


run.WSS <- function(x) {
  wss <- burden.mlogit(x, burden="WSS", ref.level=0)$p.value
  wss.pool <- burden.mlogit(x, burden="WSS", ref.level=0, group=ifelse(x@ped$pheno == 0, 0, 1))$p.value

  x.1 <- select.inds(x, x@ped$pheno != 2)
  x.2 <- select.inds(x, x@ped$pheno != 1)
  wss.cas1 <- burden.mlogit(x.1, burden="WSS", ref.level=0)$p.value
  wss.cas2 <- burden.mlogit(x.2, burden="WSS", ref.level=0)$p.value

  data.frame(wss = wss, wss.pool = wss.pool, wss.1 = wss.cas1, wss.2 = wss.cas2)
}


run.CAST <- function(x, maf.threshold) {
  cast <- burden.mlogit(x, burden="CAST", ref.level=0, maf.threshold=1)$p.value
  cast.pool <- burden.mlogit(x, burden="CAST", ref.level=0, maf.threshold=1, group=ifelse(x@ped$pheno == 0, 0, 1))$p.value

  x.1 <- select.inds(x, x@ped$pheno != 2)
  x.2 <- select.inds(x, x@ped$pheno != 1)
  cast.cas1 <- burden.mlogit(x.1, burden="CAST", ref.level=0, maf.threshold=1)$p.value
  cast.cas2 <- burden.mlogit(x.2, burden="CAST", ref.level=0, maf.threshold=1)$p.value
  
  data.frame(cast = cast, cast.pool = cast.pool, cast.cas1 = cast.cas1, cast.cas2 = cast.cas2)
}
