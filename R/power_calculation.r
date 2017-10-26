pow <- function(maf, nb.replicate, nb.controls, nb.case, prot=FALSE, filter.maf, nb.prot=0.05, nb.del=0.20){
  respow <- data.frame(test=rep(c("CAST", "WSS", "CALPHA", "BETAM", "CASTP", "WSSP", "CALPHAP", "BETAMP"), each=9), nb_gpe=rep(c(3,2), each=36), OR=rep(seq(1,5, by=0.5), 8), puissance=rep(NA, 72))
  for (OR in seq(1,5, by=0.5)){
    res <- list(CAST=rep(NA, nb.replicate), WSS=rep(NA, nb.replicate), CALPHA=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)), BETAM=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)), CASTP=rep(NA, nb.replicate), WSSP=rep(NA, nb.replicate), CALPHAP=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)), BETAMP=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)))
    for (i in 1:nb.replicate){
        maf.controls <- mafs(maf, 1)$maf.controls
    if(prot==FALSE){
        maf.case1 <- mafs(maf, sample(c(1,OR), length(maf), replace=TRUE, prob=c( (1 - nb.del), nb.del )))$maf.cases
        maf.case2 <- mafs(maf, sample(c(1,OR), length(maf), replace=TRUE, prob=c( (1 - nb.del), nb.del )))$maf.cases
    }
    else{
        maf.case1 <- mafs(maf, sample(c(1/OR,1,OR), length(maf), replace=TRUE, prob=c(nb.prot,(1 - nb.del - nb.prot), nb.del)))$maf.cases
        maf.case2 <- mafs(maf, sample(c(1/OR,1,OR), length(maf), replace=TRUE, prob=c(nb.prot,(1 - nb.del - nb.prot), nb.del)))$maf.cases
    }
      xt <- geno.simu.controls(maf.tem=maf.controls, nb.controls=nb.controls)
      xc1 <- geno.simu.case(maf.case=maf.case1, nb.case)
      xc2 <- geno.simu.case(maf.case=maf.case2, nb.case)
      xc1@ped$id <- xc1@ped$famid <- paste("C1", seq(1, nrow(xc1)), sep="_")
      xc2@ped$id <- xc2@ped$famid <- paste("C2", seq(1, nrow(xc2)), sep="_") 
      x <- rbind(xt, xc1, xc2) ; x@ped$pheno <- c(rep(0, nb.controls), rep(c(1, 2), each=nb.case)) ; x@ped$pheno <- as.factor(x@ped$pheno)
    if(filter.maf=="controls"){
        x <- select.snps(x, select.inds(x, x@ped$pheno==0)@snps$maf<0.01)
    }
    if(filter.maf=="sample"){
        x <- select.snps(x, x@snps$maf<0.01)
    }
    if(filter.maf=="any"){
        x <- select.snps(x, ( select.inds(x, x@ped$pheno==0)@snps$maf<0.01 | select.inds(x, x@ped$pheno==1)@snps$maf<0.01 | select.inds(x, x@ped$pheno==2)@snps$maf<0.01 ) )
    }
       
      res$CAST[i] <- Chi2( B_cast(x, as.factor(rep(1,ncol(x)))), x@ped$pheno)
      res$WSS[i] <- kruskal.test( WSS(x, as.factor(rep(1, ncol(x)))) , x@ped$pheno)$p.value
      res$CALPHA[i,] <- C.ALPHA.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
      res$BETAM[i,] <- Beta.M.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
      
      x@ped$pheno <- c(rep(0, nb.controls), rep(1, 2*nb.case)) ; x@ped$pheno <- as.factor(x@ped$pheno)
      res$CASTP[i] <- Chi2(B_cast(x, as.factor(rep(1,ncol(x)))), x@ped$pheno)
      res$WSSP[i] <- kruskal.test(WSS(x, as.factor(rep(1,ncol(x)))), x@ped$pheno)$p.value
      res$CALPHAP[i,] <- C.ALPHA.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
      res$BETAMP[i,] <- Beta.M.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
    }
  
  respow[which(respow$OR==OR & respow$test=="CAST"),"puissance"] <- sum(res$CAST<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="WSS"),"puissance"] <- sum(res$WSS<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="CALPHA"),"puissance"] <- sum(res$CALPHA[,"p"]<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="BETAM"),"puissance"] <- sum(res$BETAM[,"p"]<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="CASTP"),"puissance"] <- sum(res$CASTP<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="WSSP"),"puissance"] <- sum(res$WSSP<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="CALPHAP"),"puissance"] <- sum(res$CALPHAP[,"p"]<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="BETAMP"),"puissance"] <- sum(res$BETAMP[,"p"]<0.05)/nb.replicate
  }
  
  return(respow)}
  
  
  
  
  
pow2 <- function(maf, nb.replicate, nb.groups, nb.controls, nb.case, prot=FALSE, filter.maf, nb.prot=0.05, nb.del=0.20){
  respow <- data.frame(test=rep(c("CAST", "WSS", "CALPHA", "BETAM", "CASTP", "WSSP", "CALPHAP", "BETAMP"), each=9), nb_gpe=rep(c(3,2), each=36), OR=rep(seq(1,5, by=0.5), 8), puissance=rep(NA, 72))
  for (OR in seq(1,5, by=0.5)){
    res <- list(CAST=rep(NA, nb.replicate), WSS=rep(NA, nb.replicate), CALPHA=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)), BETAM=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)), CASTP=rep(NA, nb.replicate), WSSP=rep(NA, nb.replicate), CALPHAP=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)), BETAMP=data.frame(A=rep(NA, nb.replicate), B=rep(NA, nb.replicate), p=rep(NA, nb.replicate)))
    for (i in 1:nb.replicate){
        maf.controls <- mafs(maf, 1)$maf.controls
    if(prot==FALSE){
        maf.case <- replicate(mafs(maf, sample(c(1,OR), length(maf), replace=TRUE, prob=c( (1 - nb.del), nb.del )))$maf.cases, nb.groups)
    }
    else{
        maf.case <- replicate(mafs(maf, sample(c(1/OR,1,OR), length(maf), replace=TRUE, prob=c(nb.prot,(1 - nb.del - nb.prot), nb.del)))$maf.cases, nb.groups)
    }
      xt <- geno.simu.controls(maf.tem=maf.controls, nb.controls=nb.controls)
      xc <- apply(maf.case, 2, geno.simu.case(maf.case=maf.case, nb.case))
      xc@ped$id <- xc@ped$famid <- paste(rep(paste("C", seq(1, nb.groups), sep=""), nb.case), seq(1, nb.case), sep="_")
      x <- rbind(xt, xc) ; x@ped$pheno <- c(rep(0, nb.controls), rep(seq(1,nb.groups), each=nb.case)) ; x@ped$pheno <- as.factor(x@ped$pheno)
    if(filter.maf=="controls"){
        x <- select.snps(x, select.inds(x, x@ped$pheno==0)@snps$maf<0.01)
    }
    if(filter.maf=="sample"){
        x <- select.snps(x, x@snps$maf<0.01)
    }
    if(filter.maf=="any"){
        x <- select.snps(x, ( select.inds(x, x@ped$pheno==0)@snps$maf<0.01 | select.inds(x, x@ped$pheno==1)@snps$maf<0.01 | select.inds(x, x@ped$pheno==2)@snps$maf<0.01 ) )
    }
       
      res$CAST[i] <- Chi2( B_cast(x, as.factor(rep(1,ncol(x)))), x@ped$pheno)
      res$WSS[i] <- kruskal.test( WSS(x, as.factor(rep(1, ncol(x)))) , x@ped$pheno)$p.value
      res$CALPHA[i,] <- C.ALPHA.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
      res$BETAM[i,] <- Beta.M.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
      
      x@ped$pheno <- c(rep(0, nb.controls), rep(1, nb.groups*nb.case)) ; x@ped$pheno <- as.factor(x@ped$pheno)
      res$CASTP[i] <- Chi2(B_cast(x, as.factor(rep(1,ncol(x)))), x@ped$pheno)
      res$WSSP[i] <- kruskal.test(WSS(x, as.factor(rep(1,ncol(x)))), x@ped$pheno)$p.value
      res$CALPHAP[i,] <- C.ALPHA.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
      res$BETAMP[i,] <- Beta.M.2(x, centre=x@ped$pheno, region=as.factor(rep(1, ncol(x))), target=50, B.max=1000)[,c("A", "B", "p")]
    }
  
  respow[which(respow$OR==OR & respow$test=="CAST"),"puissance"] <- sum(res$CAST<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="WSS"),"puissance"] <- sum(res$WSS<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="CALPHA"),"puissance"] <- sum(res$CALPHA[,"p"]<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="BETAM"),"puissance"] <- sum(res$BETAM[,"p"]<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="CASTP"),"puissance"] <- sum(res$CASTP<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="WSSP"),"puissance"] <- sum(res$WSSP<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="CALPHAP"),"puissance"] <- sum(res$CALPHAP[,"p"]<0.05)/nb.replicate
  respow[which(respow$OR==OR & respow$test=="BETAMP"),"puissance"] <- sum(res$BETAMP[,"p"]<0.05)/nb.replicate
  }
  
  return(respow)}
  
 