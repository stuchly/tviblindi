firstleft_iter<-function(adj,iter=1000,eps=1e-6){
    tt<-as.vector(Matrix::colSums(adj))
    tt<-matrix(tt/sum(tt),nrow=1)
    j<-1
    err<-1
    while(err>eps & j<=iter){
        j<-j+1
        ttn<-tt %*% adj
        err<-sqrt(sum((tt-ttn)^2))
        tt<-ttn
    }
    print(err)
    return(as.double(tt[1,]))
}

firstleft_iter_cpp<-function(adj,iter=1000,eps=1e-6){
    tt<-as.vector(Matrix::colSums(adj))
    tt<-as.double(tt/sum(tt))
    tt<-firstleft(adj,tt,iter,eps)
    return(tt$x)
}
