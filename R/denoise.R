denoise<-function(X,IND,iter=1){
    for (i in 1:iter){
        for (j in 1:nrow(X)){
            X[j,]<-colMeans(X[IND[j,],])
        }
    }
    return(X)
}

data.core<-function(DIST,perc=0.8){
    K<-ncol(DIST)
    dens<-1/DIST[,K]
    return(which(dens>quantile(dens,prob=perc)))
}

denoiseFix<-function(X,IND,iter=1){
    XX<-X
    for (i in 1:iter){
        for (j in 1:nrow(X)){
            X[j,]<-colMeans(XX[IND[j,],])
        }
    }
    return(X)
}

magic<-function(x,iter=1,K=30,kernel="SEMer",sym="max"){
  mm<-transition.matrix(x,K=K,kernel=kernel,sym=sym,dv=0,normalize=TRUE)
  X<-x$data
  for (i in 1:iter) X<-mm%*%X
  return(as.matrix(X))
}
