.upsample<-function(data,knn,inds,N,concentration=0.1,labels=NULL){
    lab<-out<-NULL
    for (i in inds){
        ee<-matrix(.rdirichlet(N,rep(concentration,ncol(knn))),byrow = TRUE,ncol=N)
        out<-rbind(out,t(t(data[knn[i,],])%*%ee))
        if (!is.null(labels)) lab<-c(lab,rep(labels[i],N))
    }
    list(out=out,lab=lab)
}

.rdirichlet<-function (n, alpha)
{
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}


upsample<-function(x,N=5,K=30,labname="default",inds=NULL,concentration=0.1){
    if(is.null(inds)) inds<-1:nrow(x$data)
    if(is.character(inds)){
        out<-NULL
        for (j in unique(inds)){
            out<-c(out,which(x$labels[[labname]]==j))
        }
        inds<-out
        rm(out)
    }
    K<-min(K,ncol(x$KNN$INDS))
    ups<-.upsample(x$data,x$KNN$INDS[,1:K]+1,inds,N,concentration,as.character(x$labels[[labname]]))
    labels<-c(as.character(x$labels[[labname]]),ups$lab)
    data<-rbind(x$data,ups$out)
    tv<-tviblindi(data=data,labels=labels)
    return(invisible(tv))
}
