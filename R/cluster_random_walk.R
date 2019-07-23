approx_path<-function(Y,x,out=seq(0,1,by=0.01)){
    res<-matrix(NA,nrow=length(out),ncol=ncol(Y))
    for (i in 1:ncol(Y)) res[,i]<-approx(y=Y[,i],x=x,xout=out,rule=2)$y
    return(res)
}

path_dist<-function(data,path,psc,metric="max",end_pen=0){
    K<-length(path$starts)
    appsP<-list()[1:K]
    for (j in 1:K){
        if (j<K) JJ<-path$v[path$starts[j]:(path$starts[j+1]-1)] else JJ<-path$v[path$starts[j]:length(path$v)]
        appsP[[j]]<-approx_path(data[JJ,],psc[JJ])
    }
    N<-nrow(appsP[[1]])
     AA<<-appsP
    res<-matrix(NA,nrow=K,ncol=K)
    for (i in 1:(K-1)) for (j in (i+1):K){

                           ## print(approx_path(data[II,],psc[II]))
                           ## print(II)
                           ## print(psc[II])
                           ## stop()
                           appdiff<-appsP[[i]]-appsP[[j]]

                           if (metric=="euclidean")  res[i,j]<-sqrt(sum(appdiff)^2)+end_pen*sqrt(sum(appdiff[N,]^2)) else if (metric=="max") res[i,j]<-max(abs(appdiff))+end_pen*max(abs(appdiff[N,])) else if (metric=="sum") res[i,j]<-sum(abs(appdiff))+end_pen*sum(abs(appdiff[N,]))
                     }
    return(as.dist(t(res)))
}


