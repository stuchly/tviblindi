connectome<-function(x){
    if (is.null(x$metaclusters)) {
        message("Computing louvain metaclusters...\n")
        gU<-graph_from_adjacency_matrix(x$dsym,weighted=TRUE,mode="undirected")
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
    for (i in clus) g_layout<-rbind(g_layout,colMeans(x$layout[which(clus==i),]))
    aa1<-Diagonal(x=rowSums(t(agDc)+agDc)^-1)%*%agDc
    ##aa1<-Diagonal(x=rowSums(agDc)^-1)%*%agDc
    g11<-graph_from_adjacency_matrix((aa1),weighted=TRUE,mode="directed")
    E(g11)$width<-E(g11)$weight*7
    ##E(g11)$width<-E(g11)$width/max(E(g11)$width)
    E(g11)$curved=TRUE
    E(g11)$arrow.size<-0.7
    E(g11)$arrow.width<-0.8
    V(g11)$label.cex = 2
    pieD<-list()[1:length(clus)]
    for (i in clus) pieD[[i]]<-as.vector(table(x$labels[clus==i]))

    cp<-rainbow(length(unique(x$labels))+1)
    colors <- list(cp)

    colors <- rep(colors,length(clus))
    igraph::plot.igraph(g11,layout=g_layout,main="infered conectome",vertex.size=7,vertex.shape = "pie",vertex.pie=pieD,vertex.pie.color=colors)
legend("topright",legend=levels(x$labels), col=colors[[1]],pch=19,cex=1.5)

}
