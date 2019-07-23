.ps.difuse<-function(adjraw,neighbors=30,maxdim=20){
    if (neighbors>ncol(adjraw$IND[,-1])) stop ("more neighbors requested than computed!\n")
    adjraw$IND<-adjraw$IND[,1:(neighbors+1)]
    adjraw$DIST<-adjraw$DIST[,1:(neighbors+1)]
    adjraw<-knn.raw2adj(adjraw)
    spL<-sparse.Laplacian.construct(adjraw)
    .dif.res<-sparse.diffuse(spL,neigen=maxdim,maxdim=maxdim)
    return(.dif.res)
}




