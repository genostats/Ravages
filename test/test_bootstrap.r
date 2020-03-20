library(Ravages);
source("bootstrap.r")
options(width=150)
x <- as.bed.matrix(x=LCT.matrix.bed, fam=LCT.matrix.fam, bim=LCT.snps)
x@ped[,c("pop", "superpop")] <- LCT.matrix.pop1000G[,c("population", "super.population")]
x <- select.inds(x, superpop=="EUR")
x@ped$pop <- droplevels(x@ped$pop)

#Group variants within known genes
x <- set.genomic.region(x)
x1 <- filter.rare.variants(x, filter = "whole", maf.threshold = 0.025, min.nb.snps = 10)
Y <- x1@ped$pop

# premier test (pas de covariable)

SKAT(x1, Y, perm.max = 500, debug = TRUE)
SKAT.bootstrap(x1, Y, perm.max = 500, debug = TRUE)


# covariables

covar <- data.frame(sex=c(sample(0:1, sum(table(Y)[c("CEU", "GBR", "FIN")]),TRUE,c(0.2,0.8)), 
                          sample(0:1, sum(table(Y)[c("TSI", "IBS")]),TRUE,c(0.8,0.2))), 
                    u=runif(length(Y)))

Pi.matrix.LCT <- Ravages:::Pi.matrix(group=Y, data=covar, formula= ~ sex, ref.level="CEU")


x2 <- select.snps(x1, genomic.region == "LCT" & maf > 0)
x2@snps$genomic.region <- droplevels(x2@snps$genomic.region)

# on foce à en faire 5000 avec max = target = 5000
SKAT.bootstrap(x2, group=Y, Pi=Pi.matrix.LCT, X = as.matrix(covar), perm.max = 5000, perm.target = 5000, debug = TRUE)

# verification avec les variantes plus rustiques
W <- diag((1-x2@snps$maf)**24)
K <- as.matrix(x2) %*% W %*% W %*% t(as.matrix(x2))

YY <- sapply(levels(Y), function(l) as.numeric(Y == l))
Ymp <- YY - Pi.matrix.LCT
Qstat(Ymp, K, colSums(YY))

B <- bootstrap_cpp(5000, Pi.matrix.LCT, as.matrix(covar), K)
c(B$mean, B$sigma**2, B$skewness, 3 + B$kurtosis) # valeurs compatibles !!

