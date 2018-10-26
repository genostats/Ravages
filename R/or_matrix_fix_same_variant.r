OR.matrix.fix.same.variant <- function (n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro){
    if (length(OR.del) != length(OR.pro))
        stop("Dimensions mismatch")
    OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
    v <- 1:n.variants
    v.del <- sample(v, prob.del*n.variants)
    v.pro <- if(length(v.del)==0){ sample(v, prob.pro*n.variants) }else{ sample(v[-v.del], prob.pro*n.variants)}
    #On estime qu'on a toujours plus qu'un variant protecteur
    v.neutres <- if(length(v.pro)>1){ v[-c(v.del, v.pro)] }else{ v[-v.del]}
    OR.tot <- matrix(rep(NA, n.variants*nrow(OR.del)), nrow=nrow(OR.del))
    OR.tot[,v.neutres] <- 1
    if(prob.del>0){OR.tot[,v.del] <- OR[,v.del+1]}
    if(prob.pro>0){OR.tot[,v.pro] <- OR[,v.pro+1+n.variants]}
    return(OR.tot)
}
