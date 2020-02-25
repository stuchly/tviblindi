orient.sim.matrix<-function(sim,pseudotime,breaks=NULL,base=2){
    N<-length(pseudotime$res)
    oriented <- remove_back_R(sparse2triples(sim), pseudotime$res)
    if (!is.null(breaks)){
        breaks<-as.numeric(cut(as.numeric(as.factor(pseudotime$res)),breaks))
        breaks<-breaks[oriented$to]-breaks[oriented$from]
        oriented$weight<-oriented$weight/base^breaks
    }
    A<-Matrix::sparseMatrix(i=oriented$from, j=oriented$to, x=oriented$weight, dims=c(N,N))
    .DD<-Matrix::rowSums(A)
    .DD <- Matrix::Diagonal(x=.DD^-1)
    A<-.DD%*%A
}
