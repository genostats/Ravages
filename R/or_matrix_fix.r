OR.matrix.fix <- function (n.variants, n.groups, OR.del, OR.pro = 1/OR.del, prob.del, prob.pro) {
    if (length(OR.del) != length(OR.pro))
        stop("Dimensions mismatch")
    
    if(!is.matrix(OR.del)){
  	if(length(OR.del)==n.variants){
    		OR.del <- matrix(OR.del, nrow=1)
    		OR.pro <- matrix(OR.pro, nrow=1)
    	}
    	else if(length(OR.del)==n.groups){
    		OR.del <- matrix( rep_len(OR.del, n.variants*length(OR.del)), nrow=n.groups)
    		OR.pro <- matrix( rep_len(OR.pro, n.variants*length(OR.pro)), nrow=n.groups)
    	}
    	else
    		stop("OR Dimensions mismatch")
    }
    	
    
    OR <- cbind(1, OR.del, OR.pro, deparse.level = 0)
    v <- 1:n.variants
    v.del <- list()
    v.pro <- list()
    v.neutres <- list()
    #Si que des variants neutres:
    if(prob.del==0 & prob.pro==0){ 
      OR.tot <- matrix(rep(1, n.variants*nrow(OR.del)), nrow=nrow(OR.del))
    }else{
      OR.tot <- matrix(rep(NA, n.variants*nrow(OR.del)), nrow=nrow(OR.del))    
      for (i in 1:nrow(OR.del)) {
        v.del[[i]] <- sample(v, prob.del*n.variants)
        v.pro[[i]] <- if(length(v.del[[i]])==0){ sample(v, prob.pro*n.variants) }else{ sample(v[-v.del[[i]]], prob.pro*n.variants)}
        #On estime qu'on a toujours plus qu'un variant protecteur
        v.neutres[[i]] <- if(length(v.pro[[i]])>1){ v[-c(v.del[[i]], v.pro[[i]])] }else{ v[-v.del[[i]]]}
        OR.tot[i,v.neutres[[i]]] <- 1
        if(prob.del>0){OR.tot[i,v.del[[i]]] <- OR[i,v.del[[i]]+1]}
        if(prob.pro>0){OR.tot[i,v.pro[[i]]] <- OR[i,v.pro[[i]]+1+n.variants]}
      }
    }
    return(OR.tot)
}

