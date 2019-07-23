rhosigma<-function(D){
    MIN_K_DIST_SCALE<-1E-3
    ith_distM<-rowMeans(D)
    D<-D-D[,1]
    ## D<-t(t(D)/D[,ncol(D)])
    sigmaV<-C_sigma(D,log2(ncol(D)))
    ith_distM<-ith_distM*MIN_K_DIST_SCALE
    ss<-which(sigmaV<ith_distM)
    sigmaV[ss]<-ith_distM[ss]
    return(D/sigmaV)
}

symmetric_uniform_distance<-function(D){
    N<-nrow(D$IND)
    D$IND<-D$IND[,-1]
    D$DIST<-D$DIST[,-1]
    K<-ncol(D$IND)
    D$DIST<-rhosigma(D$DIST)

    adj<-sparseMatrix(i=rep(1:N,each=K),j=as.numeric(t(D$IND)+1),x=exp(-as.numeric(t(D$DIST))))
    adj<-summary(Matrix::t(adj)+adj-Matrix::t(adj)*adj)
    ## adj$x<-pmax(0,-log(adj$x)) ## solved
    adj$x<--log(adj$x)
    adj<- summary(Matrix::tril(sparseMatrix(i=adj$i,j=adj$j,x=adj$x)))
    return(sparseMatrix(i=adj$i,j=adj$j,x=adj$x))
}
