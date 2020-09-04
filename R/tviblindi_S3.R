#' Constructor of tviblindi class
#'
#' \code{tviblindi}
#' @param data a numeric matrix; the (transformed and compensated) expression data.
#' @param labels character or factor vector of the same length as \code{nrow(data)};
#' assignment of cells to population, unassigned cells should be labelled as "ungated" for better plotting.
#' Alternatively, \code{labels} can be a named list of one or label vectors.
#' @param fcs_path character (optional); path to fcs file to store the results - typically the fcs file with analysed data.
#' @param events_sel integer vector (optional); the indices of data in fcs file, defaults to 'code{1:nor(data)}.
#' @param keep_intermediate bool (default FALSE); if intermediate matrices (transition probability and spar distance matrix)
#' @param ShowAllFate bool; if all teoretical fates (see \code{Walks}) should be stored and plotted
#' should be kept during computation
#'
#' @return \code{tviblindi} returns a tviblindi class object.
#'
#' @export
tviblindi<-function(data,labels,fcs_path=NULL,events_sel=NULL,keep_intermediate=FALSE, analysis_name = paste0('tviblindi_', Sys.Date())){
    new_tviblindi(data,labels,fcs_path,events_sel,analysis_name,keep_intermediate)
}

new_tviblindi<-function(data,labels,fcs_path=NULL,events_sel=NULL,analysis_name=NULL,keep.intermediate=FALSE){
    stopifnot(is.matrix(data))
    if (!is.list(labels)) {
        labels <- list(default = labels)
    } else if (is.list(labels)) {
        stopifnot(!is.null(names(labels)))
    }
    for (idx.labels in 1:length(labels)) {
        stopifnot(length(labels[[idx.labels]])==nrow(data) && (is.factor(labels[[idx.labels]]) || is.character(labels[[idx.labels]])))
    }
    stopifnot(is.logical(keep.intermediate))
    if (is.null(events_sel)) events_sel<-1:nrow(data) ##for downsampling
    if(is.null(analysis_name)) {
        username<-Sys.info()['login']
        if(is.null(username)) username<-Sys.info()['user']
        if(is.null(username)) username<-'anonymous'
        timeanddate<-as.character(Sys.time())
        analysis_name<-paste0(username, '_', timeanddate, '_tviblindi', packageVersion('tviblindi'), '_analysis')
    }

    out<-new.env(hash=TRUE)
    out$analysis_name<-analysis_name
    out$origin<-NULL
    out$data=data
    out$denoised=NULL
    out$labels<-lapply(labels, as.factor)
    names(out$labels) <- names(labels)
    out$keep<-keep.intermediate
    out$pseudotime<-NULL
    out$filtration<-NULL
    out$boundary<-NULL
    out$reduced_boundary<-NULL
    out$walks<-NULL
    out$fates<-NULL
    out$ShowAllFates=FALSE
    out$KNN<-NULL
    out$sim<-NULL
    out$dsym<-NULL
    out$clusters<-NULL
    out$metaclusters<-NULL
    out$codes<-NULL
    out$sominfo<-NULL
    out$layout<-NULL
    out$vae=NULL
    out$events_sel=events_sel
    out$fcs=fcs_path
    structure(out,class="tviblindi")
}

print.tviblindi<-function(x){
    cat("tviblindi object\n",
        "data size: ",nrow(x$data),"\n", sep = "")
    for (idx.labels in 1:length(x$labels)) {
        cat("labels (", names(labels)[idx.labels], "): ", paste(levels(x$labels[[idx.labels]]),collapse=", "),"\n", sep = "")
    }

}

Set_origin<-function(x,...){
    UseMethod("Set_origin",x)
}

#' Sets cell-of-origin, modifies x
#'
#' \code{Set_origin}
#' @param x tviblindi class object
#' @param label character or integer; either label of population of origin (cell nearest to the mean of the population will be consider as origin)
#' or index of cell-of-origin.
#' @param labels_name name of label vector to use (only needs to be specified if multiple label vectors are present in the \code{x} and \code{label} is a string).
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
Set_origin.tviblindi<-function(x,label, labels_name = names(x$labels)[1]){
    stopifnot(length(label)==1)
    if (is.integer(label)){
        x$origin<-label
        return(invisible(x))
    }
    stopifnot(is.character(label))
    stems<-which(x$labels[[labels_name]]==label)
    x$origin <- stems[which.min(rowSums(t(t(x$data[stems, , drop = FALSE]) - colMeans(x$data[stems, , drop = FALSE]))^2))]
    return(invisible(x))
}

KNN<-function(x,...){
    UseMethod("KNN",x)
}

#' Computes KNN matrix, modifies x
#'
#' \code{KNN}
#' @param x tviblindi class object.
#' @param K integer (default K=100); number of nearest neighbors.
#' @param method character ("BT" or "annoy";default "annoy"); implements either ball tree nearest neigbor search (https://github.com/lvdmaaten/bhtsne)
#' or Approximate Nearest Neighbors Oh Yeah (https://github.com/spotify/annoy).
#' @param trees integer (default 150); number of trees for annoy - more trees more precision and more time of computation.

#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
KNN.tviblindi<-function(x,K=100,method="annoy",trees=150,denoised=FALSE){
    if (denoised & is.null(x$denoised)) stop("call Denoise() first")
    if (method=="annoy"){
        if (denoised) x$KNN<-KNN.annoy(x$denoised,K,trees) else x$KNN<-KNN.annoy(x$data,K,trees)

    } else {
        if (denoised) x$KNN<-knn.adj.raw.parallel(x$denoised, K) else x$KNN<-knn.adj.raw.parallel(x$data, K)
        x$KNN$IND<-matrix(as.integer(x$KNN$IND),ncol=ncol(x$KNN$IND))
    }
    return(invisible(x))
}

Denoise<-function(x,...){
    UseMethod("Denoise",x)
}

#' Reduce noise in data for the purpose of triangulation, modifies x
#'
#' \code{Denoise}
#' @param x tviblindi class object.
#' @param K integer (default K=30); number of neigbors to average (see details).
#' @param iter integer (default 1); number of iterations (see details).
#'
#' @details A simple noise reducution algorithm is applied - every cell coordinates are replaced by the average of \code{K} nearest neigbors.
#' This process is repeated \code{iter}-times.
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
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

#' Computes SOM clusters for triangulation, modifies x
#'
#' \code{Som}
#' @param x tviblindi class object.
#' @param xdim, ydim integer (default 25); SOM mesh size
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
Som.tviblindi<-function(x,xdim=25,ydim=25,method="som",kmeans_algorithm=c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen")){
    if (is.null(x$denoised)) {
        warning("Using original data!")
        x$denoised<-x$data
    }
    K<-xdim*ydim
    codes <-sample_points(x$denoised,K)
    ###METHOD CHANGED
    if (method!="som"){
        cl<-kmeans(x$denoised,K,centers=codes,algorithm=kmeans_algorithm)
        x$clusters<-cl$cluster
        x$codes<-cl$centers
    }
    else {
        som<-FlowSOM::SOM(x$denoised,codes=codes, xdim = xdim, ydim = ydim)
        x$clusters<-som$mapping[, 1]
        x$codes<-som$codes
        x$sominfo<-c(xdim,ydim)
    }


    return(invisible(x))
}

Filtration<-function(x,...){
    UseMethod("Filtration",x)
}

#' Computes triangulation, boundary matrix and reduced boundary matrix, modifies x
#'
#' \code{Filtration}
#' @param x tviblindi class object.
#' @param K integer (default K=30); number of nearest neighbors.
#' @param method character (only "witness" complex is implemeted); Uses Gudhi and CGAL libraries to compute witness complex
#' @param alpha double (default \code{NULL}); relaxation parameter for witness complex. If \code{NULL} mean distance to
#' K-th nearest witness is used.
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
Filtration.tviblindi<-function(x,method="witness",K=30,alpha2=NULL){
    if (!(method %in% c("traingulation","witness"))) stop("Not yet implemented")
    stopifnot(!is.null(x$codes))

    if (method=="witness"){
        xy <- FNN::get.knnx(x$codes, x$denoised, k = K)
        if (is.null(alpha2)){
            alpha2<-mean(xy$nn.dist[,K])
            cat("alpha2 = ",alpha2,"\n")
        }
        Ilist           <- split(xy$nn.index, seq(nrow(xy$nn.index)))
        Dlist           <- split(xy$nn.dist, seq(nrow(xy$nn.index)))
        x$filtration           <- witness_from_distances_cliques(Ilist, Dlist, alpha2 = alpha2, maxdimension = 1)
        x$filtration           <- create_k_skeleton(coordinates = x$codes, filtration = x$filtration, k = 2)
    }

    if (method=="traingulation"){
        ###METHOD CHANGED
        if (is.null(x$sominfo)) stop("traingulation works only with SOM clustering")
        x$filtration<-traingulation(x$codes,x$sominfo[1],x$sominfo[2])
    }

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

#' Computes pseudotime, modifies x
#'
#' \code{Pseudotime}
#' @param x tviblindi class object.
#' @param K integer (default K=30); number of nearest neighbors to compute transition matrix.
#'
#' @details Computes average distance of each cell from the cell-of-origin of all random walks in undirected graph of nearest neigbors.
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
Pseudotime.tviblindi<-function(x,K=30,nb_it=1500,iguess=NULL,eps=1e-15,kernel="Exp",kepsilon=NULL,sym="max"){
    stopifnot(!is.null(x$origin))
    if (length(x$origin)==0) stop("Origin not set!")
    if (K>dim(x$KNN$IND)[2]){
        K<-min(K,dim(x$KNN)[2])
        warning("K > dim(KNN)[2]; K<-min(K,dim(x$KNN)[2])")
    }
    d<-KofRawN(x$KNN,K)
    d  <- knn.raw2adj(d)
    dsym <- knn.spadj2sym(knn.adj2spadj(d))

    ## sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    ##METHOD CHANGED
    if (sym=="mean")
        sim <- knn.spadj.symmetrize(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    else if (sym=="prob")
        sim <- knn.spadj.symmetrize.P(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    else if (sym=="max")
        sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    else if (sym=="min"){
        d<-t(summary(dsym))
        sim <- knn.adj2spadjsim1(d, kernel = kernel,epsilon=kepsilon)
    } else stop("symmetrisation not implemented")

    x$pseudotime  <- assign_distance(sim, x$origin,weights = dsym,nb_it=nb_it,iguess=iguess,eps=eps)
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

#' Simulates random walks in directed graph, modifies x
#'
#' \code{Walks}
#' @param x tviblindi class object.
#' @param N integer (default 1000); number of walks to simulate (see details).
#' @param breaks integer (default 100); number of bins with respect to pseudotime (see details).
#' @param base double (default 1.5); penalty of jumping to far ahead in presudotime.
#' @param K integer (default 30); number of nearest neighbors to computer transition matrix.
#' @param equinumerous bool (default FALSE); simulate equal number (\code{N} for each) of walks for every potential end (see details).
#' @param to integer or character vector; indices of cells or label of target population(s) - force choice of ends (see details).
#' @param labels_name character; name of label vector to use if \code{x} has multiple label vectors and \code{to} is a string.
#' @param add bool (default FALSE); add the simulated walks to \code{x} instead of replace.
#' @param kernel character; see \code{knn.adj2spadjsim}.
#'
#' @details This method simulates random walks on directed graph (only edges pointing ahead in pseudotime are kept). To avoid short circuits
#' the pseudotime (the cells with respect to pseudotime) is divided into \code{breaks} bins and the probability of jump over k-bins ahead is penalized
#' by \code{base}^-k. Transition matrix is constructed independently of pseudotime estimation using \code{K} nearest neighbors.
#' If \code{equinumerous==TRUE} the potential ends (vertices in graph with no out-going edge) are identified and for each a subcomponent
#' of graph of vertices from which this end could be reached is used for simulation - this could be time consuming for large number of ends.
#' The same approach is used when \code{!is.null(to)}.
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
Walks.tviblindi<-function(x,N=1000,breaks=100,base=1.5,K=30, equinumerous=FALSE,to=NULL, labels_name = 'default', add=FALSE,kernel="Exp",kepsilon=NULL,sym="max"){
    if (length(x$origin)==0) stop("Origin not set!")
    add.walks<-function(x,walks){
        if(is.null(x$walks)) x$walks<-list(starts=NULL,v=NULL)
        x$walks$starts<-c(x$walks$starts,walks$starts+length(x$walks$v))
        x$walks$v<-c(x$walks$v,walks$v)
        return(invisible(0))
    }
    d<-KofRawN(x$KNN,K)
    d  <- knn.raw2adj(d)
    if (x$keep) x$dsym <- knn.spadj2sym(knn.adj2spadj(d))
    ##METHOD CHANGED
    ## sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    if (sym=="mean")
        sim <- knn.spadj.symmetrize(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    else if (sym=="prob")
        sim <- knn.spadj.symmetrize.P(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    else if (sym=="max")
        sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon))
    else if (sym=="min"){
        d<-t(summary(knn.spadj2sym(knn.adj2spadj(d))))
        sim <- knn.adj2spadjsim1(d, kernel = kernel,epsilon=kepsilon)
    } else if (sym=="none")
        sim <- knn.adj2spadjsim(d, kernel = kernel,epsilon=kepsilon)
    else stop("symmetrisation not implemented")

    if (x$keep) x$sim<-sim
    oriented.sparseMatrix <- orient.sim.matrix(sim, x$pseudotime, breaks = breaks, base = base)

    if (!equinumerous & is.null(to)){
        ## Simulate random walks
        walks           <- random_walk_adj_N_push(oriented.sparseMatrix, x$origin, N)
        if (!add){
            x$walks<-walks
        } else{
            add.walks(x,walks)
        }

    }

    fates<-NULL
    if (!is.null(to)){
        if (is.character(to)){
            to<-to[which(to %in% levels(x$labels[[labels_name]]))]
            if (length(to)==0) stop("Populations not in labels!")
            toc<-to
            to<-NULL
            for (tt in toc){
                p1<-which(x$labels[[labels_name]]==tt)
                to <- c(to,p1[which.min(rowSums(t(t(x$data[p1,]) - colMeans(x$data[p1, ]))^2))])
            }
        }
        fates <- to
        equinumerous<-TRUE
    }

    x$fates <- which(!(1:nrow(x$data) %in% Matrix::summary(oriented.sparseMatrix)$i))
    if (equinumerous){
        if (is.null(fates)) fates<-x$fates
        g<-igraph::graph_from_adjacency_matrix(oriented.sparseMatrix,weighted=TRUE,mode="directed")
        V(g)$names<-1:nrow(x$data)

        if (!add) x$walks<-NULL
        ii<-1
        for (fate in fates){
            cat("fate ",ii, " from ", length(fates),"\n")
            ii<-ii+1
            gf<-subcomponent(g,fate,"in")
            gf<-induced_subgraph(g,gf)
            nnf<-V(gf)$names
            gf<-Matrix::summary(igraph::as_adjacency_matrix(gf,attr="weight"))
            gf[,1]<-nnf[gf[,1]]
            gf[,2]<-nnf[gf[,2]]
            adjro<-sparseMatrix(i=gf$i,j=gf$j,x=gf$x,dims=dim(oriented.sparseMatrix))
            walks<-random_walk_adj_N_push(adjro, x$origin, N)
            add.walks(x,walks)
        }
    }




    return(invisible(x))
}

ToggleShowFates<-function(x,...){
    UseMethod("ToggleShowFates",x)
}

#' Helper method to set ShowAllFates, modifies x
#'
#' \code{ToggleShowFates}
#' @param x tviblindi class object.
#' @param ShowAllFates; if NULL x$ShowAllFates<-!x$ShowAllFates
#'
#' @details If \code{x$ShowAllFates==TRUE} all potential fates (see \code{Walks}) are displayed in shiny app, otherwise only ends
#' of actual simulations are displayed. \code{x$ShowAllFates==FALSE} by default.
#'
#' @return returns an invisible tviblindi class object.
#'
#' @export
ToggleShowFates.tviblindi<-function(x,ShowAllFates=NULL){
    if (is.null(ShowAllFates)) x$ShowAllFates<-!x$ShowAllFates else x$ShowAllFates<-ShowAllFates
    return(invisible(x))
}

DimRed<-function(x,...){
    UseMethod("DimRed",x)
}


##METHOD CHANGED
#' Dimensional reduction, modifies x
#'
#' \code{DimRed}
#' @param x tviblindi class object.
#' @param method character; "vaevictis" - deep (auto)encoder (combination of ideas from different papers - to be eleborated)
#' or "diffuse" - saparse diffusion maps
#' @param layout numeric matrix \code{[nrow(data),2]} (optional); if \code{!is.null(layout)} dimensional reduction is plugged-in.
#' @param dim integer (default 2); dimension of reduced data (for the moment only resonable value is 2 but the reduction should work for
#' any integer value (for "diffuse" \code{< nrow(data)})).
#' @param vsplit double (default 0.1); percentage of data used as validation step in "vaevictis".
#' @param enc_shape integer vector (default \code{c(128,128,128)}); shape (depth and wisth) of the encoder.
#' @param dec_shape integer vector (default \code{c(128,128,128)}); shape (depth and wisth) of decoder.
#' @param perplexity double (default 10.); perplexity for tsne regularisation see https://www.nature.com/articles/s41467-018-04368-5.
#' @param batch_size integer (default 512); batch size for "vaevictis" training.
#' @param epochs integer; maximum number of epochs for "vaevictis" training.
#' @param patience integer; maxim patience for for "vaevictis" training (early stopping).
#' @param ww vector double; weights for vaevictis in this order - tsne_regularisations, ivis pn loss, reconstruction error, KL divergence
#' @param margin double; ivis pn loss margin
#' @param shuffle logical; shuffle data before validation split; involves recomputation of KNN matrix
#' @param neigen integer; for "diffuse" number of eigen vectors to compute.
#' @param t double; time parameter for "diffuse", if \code{t==0} multi-time scale is used (geometric sum).
#' @param load_model character vector of 2 components; paths to files created by by x$vae$save(file1,file2) - model is loaded and applied
#' @param upsample named list \code{list(N=,takeall=)} or \code{NULL};  sample events by labels (involves recomputation of KNN matrix); affects "vaevictis" only; if NULL nothing happens, \code{N} events per label, takes all events from labels in character vector \code{takeall}
#' @param labels_name name of label vector if one is used for upsampling and \code{x} has mutliple label vectors.
#'
#' @details The pathway analysis visualisation benefits from dimensional reductions which are by definition continuous... to be elaborated
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
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
             ivis_pretrain=0,
             ww=c(10.,10.,1.,1.),
             margin=1.,
             shuffle=FALSE,
             neigen = 2,
             t = 0,
             K=30,
             load_model=NULL,
             upsample=NULL,
             labels_name = names(x$labels)[1]) {

        vae <- NULL
        if (is.null(layout)) {
            if (method[1] == "vaevictis") {
                vv = reticulate::import("vaevictis")
                if (!is.null(load_model)){
                    model <- vv$loadModel(config_file = load_model[1],weights_file = load_model[2])
                    layout <- model[[2]](x$data)
                    vae <- model
                } else {
                    if (!is.null(upsample)){
                        if (is.null(upsample$cluster)) labl<-x$labels[[labels_name]]
                        else {
                            message("running clara clustering...")
                            labl<-cluster::clara(x$data,k=upsample$cluster)$clustering
                            message("~done\n")
                        }
                        ss<-.upsample.labels(labl,N=upsample$N,takeall = upsample$takeall)
                        knn_loc<-KNN.annoy(x$data[ss,], K, 150)$IND
                        layout = vv$dimred(
                            x$data[ss,],
                            as.integer(dim),
                            vsplit,
                            enc_shape,
                            dec_shape,
                            perplexity,
                            as.integer(batch_size),
                            as.integer(epochs),
                            as.integer(patience),
                            as.integer(ivis_pretrain),
                            ww,
                            "euclidean",
                            margin,
                            K,
                            knn_loc
                        )
                    } else {
                        if (shuffle) sshuf<-sample(nrow(x$data)) else sshuf<-1:nrow(x$data)
                        if (shuffle){
                            knn.plc<-KNN.annoy(x$data[sshuf,],  K, 150)$IND
                        } else {
                            if (!is.null(x$KNN)) knn.plc<-KofRawN(x$KNN,K) else knn.plc<-KNN.annoy(x$data[sshuf,],  K, 150)$IND
                        }
                        layout = vv$dimred(
                            x$data[sshuf,],
                            as.integer(dim),
                            vsplit,
                            enc_shape,
                            dec_shape,
                            perplexity,
                            as.integer(batch_size),
                            as.integer(epochs),
                            as.integer(patience),
                            as.integer(ivis_pretrain),
                            ww,
                            "euclidean",
                            margin,
                            K,
                            knn.plc
                        )
                    }
                    x$vae <- layout[[3]]
                    #x$vae_structure<-list(config=layout[[3]]$get_config(),weights=layout[[3]]$get_weights())
                    layout <- layout[[2]](x$data)
                }


            } else if (method[1] == "diffuse") {
                if (is.null(x$KNN)) stop("Compute KNN graph first.")
                layout<-sparse.diffuse(
                    sparse.Laplacian.construct(knn.raw2adj(x$KNN)),
                    neigen = neigen,
                    t = t
                )$X

            } else if (method[1]=="umap"){
                if (!require(umap)) stop("install umap first")
                layout<-umap::umap(x$data)$layout
            }  else {
                message("Unimplemented method. Nothing done.")
            }
        }

        if (is.null(x$layout)) {
            x$layout <- list(layout)
            names(x$layout) <- paste0('1_', method[1])
        } else {
            idx.layout <- length(x$layout) + 1
            x$layout <- c(x$layout, list(layout))
            names(x$layout)[idx.layout] <- paste0(idx.layout, '_', method[1])
        }

        return(invisible(x))
    }

DownSample<-function(x,...){
    UseMethod("DownSample",x)
}

#' Downsamples data, modifies x
#'
#' \code{DownSample}
#' @param x tviblindi class object.
#' @param N integer; size of the domwnsampled data.
#' @param K integer; number of nearest neighbors to estimate density.
#' @param method character; deprecated.
#' @param e double; deprecated.
#' @param D integer; expected intrinsic dimension of the data - higher dimension underestimates the density of sparse regions.
#'
#' @details Density-based downsampling of the data, should reduce over-abundant dense populations. The labels, events_sel, data and layout are downsampled,
#' the rest is set to \code{NULL}
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
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
    x$fates<-NULL
    x$KNN<-NULL
    x$sim<-NULL
    x$dsym<-NULL
    x$clusters<-NULL
    x$metaclusters<-NULL
    x$codes<-NULL
    x$sominfo<-NULL
    x$layout<-x$layout[ss,]
    for (idx.labels in 1:length(x$labels)) {
        x$labels[[idx.labels]]<-x$labels[[idx.labels]][ss]
    }
    x$data<-x$data[ss,]
    x$events_sel<-x$events_sel[ss]
    return(invisible(x))
}

## Already generic
##METHOD CHANGED
plot.tviblindi<-function(x,pch=".",col=c("labels","pseudotime"),labels_name = names(x$labels)[1],layout_name=names(x$layout)[1],legend="bottomleft",l_cex=0.5,...){
    if (is.null(x$layout)) stop("Layout not computed!")
    if (col[1]=="pseudotime"){
        if(is.null(x$pseudotime)) stop("Pseudotime not computed!")
        psc  <- as.numeric(as.factor(x$pseudotime$res))
        psc  <- psc / max(psc)
        psc  <- psc * 10000 + 1
        col  <- greenred(10500)
        plot(x$layout[[layout_name]],col=col[psc],pch=pch,...)
    } else {
        KK<-length(levels(x$labels[[labels_name]]))
        palette(rainbow(KK))
        plot(x$layout[[layout_name]],col=x$labels[[labels_name]],pch=pch,...)
        legend(legend,legend=levels(x$labels[[labels_name]]),col=1:KK,pch=19,cex=l_cex)
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
#' Deep copy of x
#'
#' \code{Copy}
#' @param x tviblindi class object.
#'
#' @return  returns a tviblindi class object.
#'
#' @export
Copy.tviblindi<-function(x){
    y<-new.env()
    all.vars<-ls(x)
    for(var in all.vars){
        val<-x[[var]]
        y[[var]] <- if(is.environment(val)) env.copy(val,all.names) else x[[var]]
    }
    structure(y,class="tviblindi")
}

#' Computes connectome, modifies x
#'
#' \code{Pseudotime}
#' @param x tviblindi class object.
#'
#' @details See \code{connectome}
#'
#' @return  returns an invisible tviblindi class object.
#'
#' @export
Connectome<-function(x,...){
    UseMethod("Connectome",x)
}

Connectome.tviblindi<-function(x,...){
    connectome(x,...)
}







