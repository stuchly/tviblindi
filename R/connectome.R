connectome<-function(x,png="connectome.png",K=70){
    if (is.null(x$metaclusters)) {
        message("Computing louvain metaclusters...\n")
        if (is.null(x$dsym)){
            d<-KofRawN(x$KNN,K)
            d  <- knn.raw2adj(d)
            dsym <- knn.spadj2sym(knn.adj2spadj(d))
        } else dsym<-x$dsym
        gU<-graph_from_adjacency_matrix(dsym,weighted=TRUE,mode="undirected")
        x$metaclusters<-cluster_louvain(gU)$membership
        message("~Done!\n")
        rm(gU)
    }

    e.list <- cbind(x$walks$v[-c(x$walks$starts[-1]-1,length(x$walks$v))],x$walks$v[-x$walks$starts])
    e.list<-rbind(e.list,dim(x$sim))
    gD<-graph_from_edgelist(e.list,directed=TRUE)
    gD<-contract(gD,x$metaclusters)
    agDc<-summary(as_adjacency_matrix(gD))
    sel<-which(agDc$i==agDc$j)
    agDc<-agDc[-sel,]

    clus<-sort(unique(x$metaclusters))
    agDc<-sparseMatrix(i=agDc$i,j=agDc$j,x=agDc$x,dims=rep(length(clus),2))

    g_layout<-NULL
    for (i in clus) g_layout<-rbind(g_layout,colMeans(x$layout[which(x$metaclusters==i),]))
    aa1<-Diagonal(x=Matrix::rowSums(Matrix::t(agDc)+agDc)^-1)%*%agDc
    ##aa1<-Diagonal(x=rowSums(agDc)^-1)%*%agDc
    g11<-graph_from_adjacency_matrix((aa1),weighted=TRUE,mode="directed")
    E(g11)$width<-E(g11)$weight^1.2*19
    ##E(g11)$width<-E(g11)$width/max(E(g11)$width)
    E(g11)$curved=TRUE
    E(g11)$arrow.size<-2.7
    E(g11)$arrow.width<-0.7
    V(g11)$label.cex = 1
    pieD<-list()[1:length(clus)]
    for (i in clus) pieD[[i]]<-as.vector(table(x$labels[x$metaclusters==i]))

    cp<-rainbow(length(unique(x$labels))+1)
    colors <- list(cp)

    colors <- rep(colors,length(clus))
    lsize=0.4
    if (!is.null(png)){
        png(png,2000,2000)
        lsize=1.5
    }
    igraph::plot.igraph(g11,layout=g_layout,main="infered conectome",vertex.size=7,vertex.shape = "pie",vertex.pie=pieD,vertex.pie.color=colors)
    legend("topright",legend=levels(x$labels), col=colors[[1]],pch=19,cex=lsize)
     if (!is.null(png)) dev.off()
    return(invisible(x))

}
