OR.matrix.fix <- function (n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro){
    if (length(OR.del) != length(OR.pro))
        stop("Dimensions mismatch")
    OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
    v <- 1:n.variants
    v.del <- list()
    v.pro <- list()
    v.neutres <- list()
    OR.tot <- matrix(rep(NA, n.variants*nrow(OR.del)), nrow=nrow(OR.del))
    for (i in 1:nrow(OR.del)){
      v.del[[i]] <- sample(v, prob.del*n.variants)
      ifelse(length(v.del[[i]])==0, v.pro[[i]] <- sample(v, prob.pro*n.variants), v.pro[[i]] <- sample(v[-v.del[[i]]], prob.pro*n.variants))
      #On estime qu'on a toujours plus qu'un variant protecteur
      ifelse(length(v.pro[[i]])>1, v.neutres[[i]] <- v[-c(v.del[[i]], v.pro[[i]])], v.neutres[[i]] <- v[-v.del[[i]]])
      OR.tot[i,v.neutres[[i]]] <- 1
      if(prob.del>0){OR.tot[i,v.del[[i]]] <- OR[i,v.del[[i]]+1]}
      if(prob.pro>0){OR.tot[i,v.pro[[i]]] <- OR[i,v.pro[[i]]+1+n.variants]}
    }
    return(OR.tot)
    }


OR.matrix.same.fix.variant <- function (n.variants, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro){
    if (length(OR.del) != length(OR.pro))
        stop("Dimensions mismatch")
    OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
    v <- 1:n.variants
    v.del <- sample(v, prob.del*n.variants)
    ifelse(length(v.del)==0, v.pro <- sample(v, prob.pro*n.variants), v.pro <- sample(v[-v.del], prob.pro*n.variants))
    #On estime qu'on a toujours plus qu'un variant protecteur
    ifelse(length(v.pro)>1, v.neutres <- v[-c(v.del, v.pro)], v.neutres <- v[-v.del])
    OR.tot <- matrix(rep(NA, n.variants*nrow(OR.del)), nrow=nrow(OR.del))
    OR.tot[,v.neutres] <- 1
    if(prob.del>0){OR.tot[,v.del] <- OR[,v.del+1]}
    if(prob.pro>0){OR.tot[,v.pro] <- OR[,v.pro+1+n.variants]}
    return(OR.tot)
}