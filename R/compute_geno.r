geno.simu.controls <- function(maf.controls, nb.controls){
  x.controls <- matrix(rbinom(nb.controls*length(maf.controls), 2, maf.controls), byrow=TRUE, nrow=nb.controls) 
  x.controls <- as.bed.matrix(x.controls)
  x.controls@ped$pheno <- 0 ; x.controls@ped$famid <- x.controls@ped$id <- paste("T", seq(1, nb.controls), sep="")
  return(x.controls)
}


geno.simu.case <- function(maf.case, nb.case){
  x.case <- matrix(rbinom(nb.case*length(maf.case), 2, maf.case), byrow=TRUE, nrow=nb.case) 
  x.case <- as.bed.matrix(x.case)
  x.case@ped$pheno <- 1 ; x.case@ped$id <- x.case@ped$famid <- paste("C", seq(1,nb.case), sep="")
  return(x.case)
}


geno.simu <- function(maf.controls, maf.cases, nb.controls, nb.cases){
  x.controls <- matrix(rbinom(nb.controls*length(maf.controls), 2, maf.controls), byrow=TRUE, nrow=nb.controls)
  x.controls <- as.bed.matrix(x.controls)
  x.controls@ped$famid <- x.controls@ped$id <- paste("T", seq(1, nb.controls), sep="")
  
  x.cases <- matrix(rbinom(nb.cases*length(maf.cases), 2, maf.cases), byrow=TRUE, nrow=nb.cases)
  x.cases <- as.matrix(x.cases)
  x.cases@ped$famid <- x.cases@ped$id <- paste("C", seq(1, nb.cases), sep="")
  
  x <- rbind(x.controls, x.cases) ; x@ped$pheno <- rep(c(0,1), c(nb.controls, nb.cases))
  return(x)
}


random.genotypes <- function(maf, n) {
  matrix(rbinom(n*length(maf), 2, maf), byrow = TRUE, nrow = n)
}

# mafs = une matrice de mafs comme renvoyée par group.mafs
random.variants <- function(maf, n) {
  if((!is.matrix(maf) | nrow(maf) == 1) & length(n) == 1) # juste une pop
    return(as.bed.matrix(random.genotypes(maf, n)))

}
