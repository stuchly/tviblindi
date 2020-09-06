.upsample<-function(data,knn,inds,N){
    out<-NULL
    for (i in 1:length(inds)){
        ee<-matrix(runif(N*ncol(knn)),ncol=N)
        ee<-ee/colSums(ee)
        out<-rbind(out,data[knn[i,],]%*%ee)
    }
    out
}
