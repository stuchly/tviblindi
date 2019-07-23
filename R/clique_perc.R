clique_comunity<-function(graph,k){
    clq<-cliques(graph,min=k,max=k)
    N<-length(clq)
    ij<-connect_cliques(as.matrix(do.call("rbind",clq)))
    ij<-sparseMatrix(i=ij[[1]]+1,j=ij[[2]]+1,x=1,dims=c(N,N))
    clq.graph<-simplify(graph_from_adjacency_matrix(ij,mode="max"))

    V(clq.graph)$name<-seq_len(vcount(clq.graph))

    comps<-decompose.graph(clq.graph)

    lapply(comps,function(x) {unique(unlist(clq[V(x)$name]))})
}


clique_comunity_R<-function(graph,k){
    clq<-cliques(graph,min=k,max=k)
    edges<-c()
    N<-length(clq)
    print(N)
    ## e <- sparseMatrix(dims = c(N,N), i={}, j={},x=0)
    for (i in 1:(N-1)){
        print(i)
        for (j in (i+1):N){
            if (length(intersect(clq[[i]],clq[[j]]))==(k-1)) edges<-c(edges,c(i,j))
        }
    }
    ## return(e)
    clq.graph<-simplify(graph(edges))

    V(clq.graph)$name<-seq_len(vcount(clq.graph))

    comps<-decompose.graph(clq.graph)

    lapply(comps,function(x) {unique(unlist(clq[V(x)$name]))})
}
