traingulationTDA<-function(clusters,xdim,ydim) create_k_skeleton(TDA::alphaComplexFiltration(expand.grid(seq_len(xdim), seq_len(ydim))),clusters,2,type="alpha")

traingulation<-function(clusters,xdim,ydim,method="manhattan",k=10){
    grid<-expand.grid(seq_len(xdim), seq_len(ydim))
    grid<-as.matrix(dist(grid,method=method))
    th<-sort(grid[,1])[k+1]
    grid[grid>th]<-0
    ## grid<-Matrix(grid,sparse=TRUE)
    grid<-as(grid,"dgCMatrix")
    ## print(class(grid))
    ## return(grid)
    flt<-rips_from_spadj(grid, maxdimension = 1, nbVertices = dim(grid)[1], threshold = max(grid))
    create_k_skeleton(flt,clusters,2,type="alpha")
}
