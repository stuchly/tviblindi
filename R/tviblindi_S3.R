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
    out$vae_predict=NULL
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

Walks.tviblindi<-function(x,N=1000,breaks=100,base=1.5,K=30){
    
    d<-KofRawN(x$KNN,K)
    d  <- knn.raw2adj(d)
    x$dsym <- knn.spadj2sym(knn.adj2spadj(d))
    x$sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = "Exp"))
    oriented.sparseMatrix <- orient.sim.matrix(x$sim, x$pseudotime, breaks = breaks, base = base)

    ## Simulate random walks
    x$walks           <- random_walk_adj_N_push(oriented.sparseMatrix, x$origin, N)
    return(invisible(x))
}

DimRed<-function(x,...){
    UseMethod("DimRed",x)
}

DimRed.tviblindi <-
    function(x,
             method = c("vaevictis", "diffuse"),
             layout = NULL,
             dim = 2,
             vsplit = 0.1,
             enc_shape = c(128, 128, 128),
             dec_shape = c(128, 128, 128),
             perplexity = 10.,
             batch_size = 512L,
             epochs = 100L,
             patience = 0L,
             alpha = 10.,
             neigen = 2,
             t = 0) {
        if (!is.null(layout)) {
            x$layout <- layout
            return(invisible(x))
        }
        if (method[1] == "vaevictis") {
            vv = reticulate::import("vaevictis")
            layout = vv$dimred(
                x$data,
                as.integer(dim),
                vsplit,
                enc_shape,
                dec_shape,
                perplexity,
                as.integer(batch_size),
                as.integer(epochs),
                as.integer(patience),
                alpha
            )
            x$vae_predict = layout[[2]]
            x$layout <- layout[[1]]
        } else if (method[1] == "diffuse") {
            if (is.nul(x$KNN)) stop("Compute KNN graph first.")
            x$layout<-sparse.diffuse(
                sparse.Laplacian.construct(knn.raw2adj(x$KNN)),
                neigen = neigen,
                t = t
            )$X
            
        } else {
            message("Unimplemented method. Nothing done.")
        }
        return(invisible(x))
    }

## Already generic
plot.tviblindi<-function(x,pch=".",col=c("labels","pseudotime")){
    if (is.null(x$layout)) stop("Layout not computed!")
    if (col=="pseudotime"){
        if(is.null(x$pseudotime)) stop("Pseudotime not computed!")
        psc  <- as.numeric(as.factor(x$pseudotime$res))
        psc  <- psc / max(psc)
        psc  <- psc * 10000 + 1
        col  <- greenred(10500)
        plot(x$layout,col=col[psc],pch=pch)
    } else {
        plot(x$layout,col=x$labels,pch=pch)
    }
}
mock_pass<-function(x){
    nx<-deparse(substitute(x))
    saveRDS(nx,"nx")
    nz<-readRDS("nx")
    z<-get(nz)
    return(z)
}









