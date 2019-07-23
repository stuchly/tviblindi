sample_adjact<-function(g,origin,samplesize=500){
    N<-gorder(g)
    waypoints<-as.character(c(origin,sample(x=(1:N)[-origin],size=samplesize-1)))
    print(waypoints)
    adjact<-matrix(0,samplesize,samplesize)
    for (i in 1:(samplesize-1)){
        print(i)
        for (j in (i+1):samplesize) {
            if (length(shortest_paths(from=waypoints[i],to=waypoints[j],graph=g,output = "vpath")$vpath[[1]])>0 | length(shortest_paths(from=waypoints[j],to=waypoints[i],graph=g,output="vpath")$vpath[[1]])>0) adjact[i,j]<-adjact[j,i]<-1


        }
    }
    return(list(g=graph_from_adjacency_matrix(adjact,mode="undirected"),waypoints=waypoints))
}
