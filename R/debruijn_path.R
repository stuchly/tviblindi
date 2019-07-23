kmers<-function (x, k){
    K<-length(x)-k+1
    return(matrix(x[rep(1:K,each=k)+rep(1:k,K)-1],ncol=k,byrow=TRUE))
}
## klmers<-function (x, k,l){
##     L<-seq(1,length(x)-k+1,by=l)
##     print(L)
##     return(rep(L,each=k)+rep(1:k,length(L))-1)
##     return(matrix(x[rep(L,each=k)+rep(1:k,length(L))-1],ncol=k,byrow=TRUE))
## }
path_kmers<-function(paths,k){
    ind<-paths$starts
    paths<-paths$v
    res<-NULL
    for (i in 1:(length(ind)-1)){
        res<-rbind(res,kmers(paths[ind[i]:ind[i+1]],k))
    }
    res<-rbind(res,kmers(paths[ind[i]:length(paths)],k))
    return(res)
}

kmers_identical<-function(kmerA,kmerB,graph,order=1){
    if (identical(kmerA,kmerB)) return(TRUE)
    ##if (all(kmerA %in% unique(unlist(ego(graph,order,kmerB))))) return(TRUE)
    if (length(intersect(kmerA,(unlist(ego(graph,order,kmerB)))))==length(kmerA)) return(TRUE)
    ##if (length(setdiff(kmerA,unlist(ego(graph,order,kmerB))))==0) return(TRUE)
    FALSE

}

kmers_identical2<-function(kmerA,kmerB,adj){
    if (identical(kmerA,kmerB)) return(TRUE)
    if (length(intersect(kmerA,as.vector(adj[kmerB,])))==length(kmerA) || length(intersect(kmerB,as.vector(adj[kmerA,])))==length(kmerB)) return(TRUE)
    FALSE

}

kmers_identical3<-function(kmerA,kmerB,adj){
    if (identical(kmerA,kmerB)) return(TRUE)
    neigb<-union(intersect(kmerA,as.vector(adj[kmerB,])),intersect(kmerB,as.vector(adj[kmerA,])))
    if (length(intersect(kmerA,neigb))==length(kmerA)) return(TRUE)
    FALSE

}

build_graph<-function(kmers,adj){
    nodes<-1
    for (i in 2:nrow(kmers)){
        if (i %% 100 ==0 ) print(paste(i,length(nodes),sep=":"))
        found<-FALSE
        for (nn in nodes){
            if (kmers_identical2(kmers[i,],kmers[nn,],adj)) {found<-TRUE;break}
        }
        if (!found) nodes<-c(nodes,i)
    }
    return(nodes)
}
