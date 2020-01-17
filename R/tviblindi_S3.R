tviblindi<-function(data,labels){
    new_tviblindi(data,labels)
}

new_tviblindi<-function(data,labels,keep.intermediate=FALSE){
    stopifnot(is.matrix(data))
    stopifnot(length(labels)==nrow(data) && (is.factor(labels) || is.character(labels)))
    stopifnot(is.logical(keep.intermediate))

    out<-new.env(hash=TRUE)
    out$origin<-NULL
    out$data=data
    out$labels<-as.factor(labels)
    out$keep<-keep.intermediate
    out$pseudotime<-NULL
    out$filtration<-NULL
    out$boundary<-NULL
    out$reduced_boundary<-NULL
    out$walks<-NULL
    out$KNN<-NULL
    out$sim<-NULL
    out$dsim<-NULL
    out$clusters<-NULL
    out$codes<-NULL
    out$layout<-NULL
    structure(out,class="tviblindi")
}

Set_origin<-function(x,...){
    UseMethod("Set_origin",x)
}

Set_origin.tviblindi<-function(x,label){
    stopifnot(length(label)==1)
    if (is.integer(label)){
        x$origin<-label
        return(invisible(x))
    }
    stopifnot(is.character(label))
    stems<-which(x$labels==label)
    x$origin <- stems[which.min(rowSums(t(t(x$data[stems,]) - colMeans(x$data[stems, ]))^2))]
    return(invisible(x))
}

KNN<-function(x,...){
    UseMethod("KNN",x)
}

KNN.tviblindi<-function(x,K=100){
    x$KNN<-knn.adj.raw.parallel(x$data, K)
    return(invisible(x))
}

Denoise<-function(x,...){
    UseMethod("Denoise",x)
}

Denoise.tviblindi<-function(x,K=30,iter=1){
    stopifnot(!is.null(x$KNN))

    if (K>dim(x$KNN$IND)[2]){
        K<-min(K,dim(x$KNN)[2])
        warning("K > dim(KNN)[2]; K<-min(K,dim(x$KNN)[2])")
    }
        x$denoised<-denoise(x$data,x$KNN$IND[,2:K]+1,iter=iter)
    return(invisible(x))
}

Som<-function(x,...){
    UseMethod("Som",x)
}

Som.tviblindi<-function(x,xdim=25,ydim=25){
    if (is.null(x$denoised)) {
        warning("Using original data!")
        x$denoised<-x$data
    }
    K<-xdim*ydim
    codes <-sample_points(x$denoised,K)
    som<-FlowSOM::SOM(x$denoised,codes=codes, xdim = xdim, ydim = ydim)
    x$clusters<-som$mapping[, 1]
    x$codes<-som$codes
    return(invisible(x))
}

Filtration<-function(x,...){
    UseMethod("Filtration",x)
}

Filtration.tviblindi<-function(x,method="witness",K=30,alpha2=10){
    if (method!="witness") stop("Not yet implemented")
    stopifnot(!is.null(x$codes))
    xy <- FNN::get.knnx(x$codes, x$denoised, k = K)
    Ilist           <- split(xy$nn.index, seq(nrow(xy$nn.index)))
    Dlist           <- split(xy$nn.dist, seq(nrow(xy$nn.index)))
    x$filtration           <- witness_from_distances_cliques(Ilist, Dlist, alpha2 = alpha2, maxdimension = 1)
    x$filtration           <- create_k_skeleton(coordinates = x$codes, filtration = x$filtration, k = 2)
    x$boundary<-build_boundaryC(x$filtration)
    x$reduced_boundary<-reduce_boundary(x$boundary)
    return(invisible(x))
}

## Boundary<-function(x,...){
##     UseMethod("Boundary",x)
## }

## Boundary.tviblindi<-function(x){
##     stopifnot(!is.null(x$filtration))
##     x$boundary<-build_boundaryC(x$filtration)
##     x$reduced_boundary<-reduce_boundary(x$boundary)
##     return(invisible(x))
## }

Pseudotime<-function(x,...){
    UseMethod("Pseudotime",x)
}

Pseudotime.tviblindi<-function(x,K=30){
    stopifnot(!is.null(x$origin))
    if (K>dim(x$KNN$IND)[2]){
        K<-min(K,dim(x$KNN)[2])
        warning("K > dim(KNN)[2]; K<-min(K,dim(x$KNN)[2])")
    }
    d<-KofRawN(x$KNN,K)
    d  <- knn.raw2adj(d)
    x$dsym <- knn.spadj2sym(knn.adj2spadj(d))
    x$sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = "Exp"))
    x$pseudotime  <- assign_distance(x$sim, x$origin,weights = x$dsym)
    cat("Pseudotime error:", x$pseudotime$error, "\n")
    return(invisible(x))
}


Walks<-function(x,...){
    UseMethod("Walks",x)
}

Walks.tviblindi<-function(x,N=1000,breaks=100,base=1.5){

    oriented.sparseMatrix <- orient.sim.matrix(x$sim, x$pseudotime, breaks = breaks, base = base)

    ## Simulate random walks
    x$walks           <- random_walk_adj_N_push(oriented.sparseMatrix, x$origin, N)
    return(invisible(x))
}

DimRed<-function(x,...){
    UseMethod("DimRed",x)
}

DimRed.tviblindi<-function(x,layout){##for consistency only for now
    tv1$layout<-layout
    return(invisible(x))
}

mock_pass<-function(x){
    nx<-deparse(substitute(x))
    saveRDS(nx,"nx")
    nz<-readRDS("nx")
    z<-get(nz)
    return(z)
}









