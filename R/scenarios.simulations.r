#SVDR: two different groups of cases, needs to different prevalences
#otherwise: same thresholds
#SVDR <- function(R, size = c(1000,500,500), replicates = 10) {

#  s1 <- R$thresholds[1]
#  s2 <- R$thresholds[2]

#  x <- Ravages:::random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, s1, s2), c(Inf, s2, Inf), size, replicates)
#}

#DG: same group as the controls
DG <- function(R, size = c(1000, 500, 500), replicates = 10) {

  s <- R$thresholds[1]

  x <- random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, -Inf, s), c(Inf, Inf, Inf), size, replicates)
}

#SVSR: two groups of cases identical
SVSR <- function(R, size = c(1000, 500, 500), replicates = 10) {

  s <- rep(R$thresholds[1], length(size)-1)

  x <- random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, s), rep(Inf, length(size)), size, replicates)
}


#DVDR: needs two lists of parameters
#DVDR <- function(R1, R2, size = c(1000, 500, 500), replicates = 10) {

#  s1 <- R1$thresholds[1]
#  x1 <- Ravages:::random.bed.matrix.haplos.thresholds(R1$haplos, R1$burdens, sd = sqrt(1-R1$h2), c(-Inf, s1), c(Inf, Inf), size[1:2], replicates)

#  s2 <- R2$thresholds[1]
#  x2 <- Ravages:::random.bed.matrix.haplos.thresholds(R2$haplos, R2$burdens, sd = sqrt(1-R2$h2), c(-Inf, s2), c(Inf, Inf), c(0,size[3]), replicates)
#  x2@ped$pheno <- 2 
#  x <- rbind(x1,x2)
#}

DVDR <- function(R, size = c(1000, 500, 500), replicates = 10) {

  s <- as.numeric(lapply(R, function(z) z$threshold[1]))
  h2 <- as.numeric(lapply(R, function(z) z$h2))

  x1 <- random.bed.matrix.haplos.thresholds(R[[1]]$haplos, R[[1]]$burdens, sd = sqrt(1-h2[1]), c(-Inf, s[1]), c(Inf, Inf), size[1:2], replicates)
  x.temp <- lapply(2:length(s), function(z) random.bed.matrix.haplos.thresholds(R[[z]]$haplos, R[[z]]$burdens, sd = sqrt(1-h2[z]), c(-Inf, s[z]), c(Inf, Inf), c(0, size[z+1]), replicates))

  x <- do.call(rbind, c(x1,x.temp))

  x@ped$pheno <- rep(1:length(size), size)
  return(x)
}


#DVSR: different variants but same risks
DVSR <- function(R, size = c(1000, 500, 500), replicates = 10) {
  s <- R$thresholds[1]
  x1 <- random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, s), c(Inf, Inf), size[1:2], replicates)
  x.temp <- lapply(2:(length(size)-1), function(z) random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, s), c(Inf, Inf), c(0, size[z+1]), replicates))
  #x2 <- Ravages:::random.bed.matrix.haplos.thresholds(R$haplos, R$burdens, sd = sqrt(1-R$h2), c(-Inf, s), c(Inf, Inf), c(0,size[3]), replicates)
  #x2@ped$pheno <- 2 
  #x <- rbind(x1,x2)

  x <- do.call(rbind, c(x1,x.temp))
  x@ped$pheno <- rep(1:length(size), size)
  return(x)
}

