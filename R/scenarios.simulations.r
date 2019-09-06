#SVDR: two different groups of cases, needs to different prevalences
#otherwise: same thresholds
SVDR <- function(R, size = c(1000,500,500), replicates = 10) {

  s1 <- R$thresholds[1]
  s2 <- R$thresholds[2]

  x <- Ravages:::random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, s1, s2), c(Inf, s2, Inf), size, replicates)
}

#DG: same group as the controls
DG <- function(R, size = c(1000, 500, 500), replicates = 10) {

  s <- R$thresholds[1]

  x <- Ravages:::random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, -Inf, s), c(Inf, Inf, Inf), size, replicates)
}

#SVSR: two groups of cases identical
SVSR <- function(R, size = c(1000, 500, 500), replicates = 10) {

  s <- R$thresholds[1]

  x <- Ravages:::random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, s, s), c(Inf, Inf, Inf), size, replicates)
}


#DVSG: needs two lists of parameters
DVSG <- function(R1, R2, size = c(1000, 500, 500), replicates = 10) {

  s1 <- R1$thresholds[1]
  x1 <- Ravages:::random.bed.matrix.haplos.thresholds(R1$haplos, R1$burdens, sd = sqrt(1-R1$h2), c(-Inf, s1), c(Inf, Inf), size[1:2], replicates)

  s2 <- R2$thresholds[1]
  x2 <- Ravages:::random.bed.matrix.haplos.thresholds(R2$haplos, R2$burdens, sd = sqrt(1-R2$h2), c(-Inf, s2), c(Inf, Inf), c(0,size[3]), replicates)
  x2@ped$pheno <- 2 
  x <- rbind(x1,x2)
}

