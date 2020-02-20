tviblindi<-function(data,labels,fcs_path=NULL,events_sel=NULL,keep_intermediate=FALSE){
    new_tviblindi(data,labels,fcs_path,events_sel,keep_intermediate)
}

new_tviblindi<-function(data,labels,fcs_path=NULL,events_sel=NULL,keep.intermediate=FALSE){
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
    out$dsym<-NULL
    out$clusters<-NULL
    out$metaclusters<-NULL
    out$codes<-NULL
    out$layout<-NULL
    out$vae_predict=NULL
    out$events_sel=events_sel
    out$fcs=fcs_path
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

KNN.tviblindi<-function(x,K=100,method="BT",trees=150){
  if (method=="annoy"){
    x$KNN<-KNN.annoy(x$data,K,trees)
  } else {
    x$KNN<-knn.adj.raw.parallel(x$data, K)
  }
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
    dsym <- knn.spadj2sym(knn.adj2spadj(d))
    sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = "Exp"))
    x$pseudotime  <- assign_distance(sim, x$origin,weights = dsym)
    cat("Pseudotime error:", x$pseudotime$error, "\n")
    if (x$keep) {
        x$sim<-sim
        x$dsym<-dsym
    }
    return(invisible(x))
}


Walks<-function(x,...){
    UseMethod("Walks",x)
}

Walks.tviblindi<-function(x,N=1000,breaks=100,base=1.5,K=30){

    d<-KofRawN(x$KNN,K)
    d  <- knn.raw2adj(d)
    if (x$keep) x$dsym <- knn.spadj2sym(knn.adj2spadj(d))
    sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = "Exp"))
    oriented.sparseMatrix <- orient.sim.matrix(sim, x$pseudotime, breaks = breaks, base = base)
    if (x$keep) x$sim<-sim
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
            if (is.null(x$KNN)) stop("Compute KNN graph first.")
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

DownSample<-function(x,...){
    UseMethod("DownSample",x)
}
DownSample.tviblindi<-function(x,N=10000,K=10,method="default",e=1.,D=2){
    if (is.null(x$KNN)) stop("Compute KNN first.")
    N=min(nrow(x$data),N)
    if (method=="naive")
        ss<-sample(x=1:nrow(x$data),prob=x$KNN$DIST[,K]^e,size=N,replace=FALSE)
    else
        if (method=="exp") ss<-sample(x=1:nrow(x$data),prob=exp(x$KNN$DIST[,K]^e),size=N,replace=FALSE)
    else {
        k<-ncol(x$data)
        Vd<-pi^(D/2)/gamma(D/2+1)
        dens<-K/(nrow(x$data)*Vd)
        dens<-dens/x$KNN$DIST[,K]^D
        dens<-dens^e
        density_s <- sort(dens)
		cdf<- rev(cumsum(1.0/rev(density_s)))
        boundary <- N/cdf[1]

		if (boundary > density_s[1]) {  # Boundary actually falls amongst densities present
			targets <- (N-1:length(density_s)) / cdf
			boundary <- targets[which.min(targets-density_s > 0)]
		}
		ss<-which(boundary/dens > runif(length(dens)))
    }

    x$pseudotime<-NULL
    x$filtration<-NULL
    x$boundary<-NULL
    x$reduced_boundary<-NULL
    x$walks<-NULL
    x$KNN<-NULL
    x$sim<-NULL
    x$dsym<-NULL
    x$clusters<-NULL
    x$metaclusters<-NULL
    x$codes<-NULL
    x$layout<-x$layout[ss,]
    x$labels<-x$labels[ss]
    x$data<-x$data[ss,]
    return(invisible(x))
}

## Already generic
plot.tviblindi<-function(x,pch=".",col=c("labels","pseudotime"),legend="bottomleft",l_cex=0.5){
    if (is.null(x$layout)) stop("Layout not computed!")
    if (col[1]=="pseudotime"){
        if(is.null(x$pseudotime)) stop("Pseudotime not computed!")
        psc  <- as.numeric(as.factor(x$pseudotime$res))
        psc  <- psc / max(psc)
        psc  <- psc * 10000 + 1
        col  <- greenred(10500)
        plot(x$layout,col=col[psc],pch=pch)
    } else {
        KK<-length(levels(x$labels))
        palette(rainbow(KK))
        plot(x$layout,col=x$labels,pch=pch)
        legend(legend,legend=levels(x$labels),col=1:KK,pch=19,cex=l_cex)
        palette("default")
    }
}
mock_pass<-function(x){
    nx<-deparse(substitute(x))
    saveRDS(nx,"nx")
    nz<-readRDS("nx")
    z<-get(nz)
    return(z)
}

Copy<-function(x,...){
    UseMethod("Copy",x)
}


##Adapted for Rfast package
Copy.tviblindi<-function(x){
    y<-new.env()
    all.vars<-ls(x)
    for(var in all.vars){
        val<-x[[var]]
        y[[var]] <- if(is.environment(val)) env.copy(val,all.names) else x[[var]]
    }
    structure(y,class="tviblindi")
}

Connectome<-function(x,...){
    UseMethod("Connectome",x)
}

Connectome.tviblindi<-function(x,...){
    connectome(x)
}







