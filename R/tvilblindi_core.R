## knn.construct<-function(XX,K,metric="euclid"){
##     D<-ncol(XX)
##     N<-nrow(XX)
##     print(N)
##     print(D)
##     res<-NULL
##     if (metric=="euclid"){ res<-.Call("_tviblindi_createKNNgraph_eu",X=as.numeric(t(XX)),D=D,K=K+1,N=N)}
##     else{
##         if (metric=="cosine") res<-.Call("_tviblindi_createKNNgraph_cos",X=as.numeric(t(XX)),D=D,K=K+1,N=N)
##     }

##     adj<-rbind(rep(1:N,each=K),as.numeric(t(res$IND[,-1])),as.numeric(t(res$DIST[,-1])))

##     gg<-graph_from_adjacency_matrix(sparseMatrix(i=adj[1,],j=adj[2,]+1,x=adj[3,],dims=c(N,N)),weighted=TRUE)
##     return(list(knn=res,knn.adj=adj,knn.graph=gg))

## }

## lofknn.construct<-function(knn,l,N){
##     IND<-knn$knn$IND[,-1]
##     DIST<-knn$knn$DIST[,-1]
##     .ind<-.dist<-matrix(NA,nrow(IND),l*N)
##     ii<-rep(1:nrow(IND),each=l)
##     for (i in 1:nrow(IND)){
##         if (i %% 10000==0) print(i)
##         Lind<-as.numeric(replicate(n=N,sample(x=1:ncol(IND),size=l)))
##         .ind[i,]<-IND[i,Lind]
##         .dist[i,]<-DIST[i,Lind]
##     }
##     gg<-list()[1:N]

##     for (i in 1:N){
##         gg[[i]]<-graph_from_adjacency_matrix(sparseMatrix(i=ii,j=as.numeric(t(.ind[,(l*(i-1)+1):(i*l)]+1)),x=as.numeric(t(.dist[,(l*(i-1)+1):(i*l)])),dims=c(nrow(IND),nrow(IND))),weighted=TRUE)

##     }

##     return(gg)
## }

#' Create a sparse adjacency matrix object
#'
#' \code{knn.adj.construct} generates an adjacency matrix object for a K-nearest-neighbours graph. Nearest neighbour indices and neighbour distances are found. Nearest neighbour search is implemented via a vantage-point tree algorithm (used in Van der Maaten, 2014).
#'
#' @param XX matrix of coordinates. Rows correspond to vertices.
#' @param K number of nearest neighbours.
#' @param metric distance metric. \code{euclid} or \code{cosine}.
#'
#' @return \code{knn.adj.construct} returns a matrix consisting of 3 rows. 1st row contains (1-based) vertex indices. 2nd row contains indices of nearest neighbours. 3rd row contains distances between data points in rows 1 and 2.
#'
#' @references
#' Van der Maaten, L (2014). Accelerating t-SNE using Tree-Based Algorithms. Journal of Machine Learning Research 15 (2014) 1-21
#'
#' @export
knn.adj.construct<-function(XX,K,metric="euclid"){
    D<-ncol(XX)
    N<-nrow(XX)
    print(N)
    print(D)
    res<-NULL
    if (metric=="euclid"){ res<-.Call("_tviblindi_createKNNgraph_eu",X=as.numeric(t(XX)),D=D,K=K+1,N=N)}
    else{
        if (metric=="cosine") res<-.Call("_tviblindi_createKNNgraph_cos",X=as.numeric(t(XX)),D=D,K=K+1,N=N)
    }

    adj<-rbind(rep(1:N,each=K),as.numeric(t(res$IND[,-1])+1),as.numeric(t(res$DIST[,-1])))
    return(adj)
}

#' Create an adjacency matrix
#'
#' \code{knn.adj.raw} generates an adjacency matrix for a K-nearest-neighbours graph. A matrix of neighbour indices and a matrix of neighbour distances are created. Nearest neighbour search is implemented via a vantage-point tree algorithm (Van der Maaten, 2013).
#'
#' @param XX matrix of coordinates. Rows correspond to vertices.
#' @param K number of nearest neighbours.
#' @param metric distance metric. \code{euclid} or \code{cosine}.
#'
#' @return
#' \code{knn.adj.raw} returns a list of matrix \code{IND} and matrix \code{DIST}.
#'
#' \code{IND} contains 0-based indices of each vertex (first column) and its \code{K} nearest neighbours, in the following columns.
#'
#' \code{DIST} contains the corresponding distances.
#'
#' @references
#'
#'
#' @examples
#' # create an artificial data set of two clusters in 2-dimensional space
#' coord <- cbind(c(rnorm(1000,100,20), rnorm(500, 20, 5)), c(rnorm(1000,100,20), rnorm(500, 120, 5)))
#' coord.adj <- knn.adj.raw(coord, K=5)
#'
#' # pick some two vertices in the kNN graph and show them and their neighbours on a plot
#' plot(coord)
#' points(coord[coord.adj$IND[1,]+1,], col="blue", pch=20)
#' points(coord[coord.adj$IND[2,]+1,], col="red", pch=20)
#'
#' @export
knn.adj.raw<-function(XX,K,metric="euclid"){
    D<-ncol(XX)
    N<-nrow(XX)
    print(N)
    print(D)
    res<-NULL
    if (metric=="euclid"){ res<-.Call("_tviblindi_createKNNgraph_eu",X=as.numeric(t(XX)),D=D,K=K+1,N=N)}
    else{
        if (metric=="cosine") res<-.Call("_tviblindi_createKNNgraph_cos",X=as.numeric(t(XX)),D=D,K=K+1,N=N)
    }

    return(res)
}

raw2cknn<-function(adjraw,K=ncol(adjraw$IND)-1){
    if (K>ncol(adjraw$DIST)) stop("K> ncol(adjraw)")

    KK<-sqrt(adjraw$DIST[,K])
    KK1<-rep(KK,each=ncol(adjraw$DIST))
    adjraw$DIST<-t(t(adjraw$DIST)/KK1)
    rm(KK1)
    adjraw$DIST<-t(t(adjraw$DIST)/KK[as.integer(t(adjraw$IND)+1)])
    return(adjraw)

}

reorder.raw<-function(adjraw){
    N<-nrow(adjraw$IND)
    for (i in 1:N){
        ord<-order(adjraw$DIST[i,])
        adjraw$DIST[i,]<-adjraw$DIST[i,ord]
        adjraw$IND[i,]<-adjraw$IND[i,ord]
    }
    return(adjraw)
}

knn.raw2jacc<-function(adjraw){
    adjraw<-jaccard_coeff(adjraw$IND[,-1]+1)
    adjraw<-adjraw[adjraw[,3]>0,]
    adjraw <- as.data.frame(adjraw)
    colnames(adjraw)<- c("from","to","weight")
    g <- igraph::graph.data.frame(adjraw, directed=FALSE)
    return(as_adjacency_matrix(g))
}

#' Convert adjacency matrix to a sparseMatrix object
#'
#' \code{knn.raw2adj} converts an adjacency matrix, as generated by \code{knn.adj.raw} to a \code{sparseMatrix} object, as generated by \code{knn.adj.construct}.
#'
#' @param adjraw adjacency matrix.
#'
#' @return \code{knn.raw2adj} returns a matrix consisting of 3 rows. 1st row contains (1-based) vertex indices. 2nd row contains indices of nearest neighbours. 3rd row contains distances between data points in rows 1 and 2.
#'
#' @examples
#' # create an artificial data set of two clusters in 2-dimensional space
#' coord <- cbind(c(rnorm(1000,100,20), rnorm(500, 20, 5)), c(rnorm(1000,100,20), rnorm(500, 120, 5)))
#'
#' coord.adj <- knn.adj.raw(coord, K=5)
#' coord.adj.sparse <- knn.raw2adj(coord.adj)
#'
#' @export
knn.raw2adj<-function(adj){
    N<-nrow(adj$IND)
    K<-ncol(adj$IND)-1
    adj<-rbind(rep(1:N,each=K),as.numeric(t(adj$IND[,-1])+1),as.numeric(t(adj$DIST[,-1])))
}

knn.spadj.construct<-function(XX,K,metric="euclid"){
    D<-ncol(XX)
    N<-nrow(XX)

    print(N)
    print(D)
    res<-NULL
    if (metric=="euclid"){ res<-.Call("_tviblindi_createKNNgraph_eu",X=as.numeric(t(XX)),D=D,K=K+1,N=N)}
    else{
        if (metric=="cosine") res<-.Call("_tviblindi_createKNNgraph_cos",X=as.numeric(t(XX)),D=D,K=K+1,N=N)
    }

    adj<-rbind(rep(1:N,each=K),as.numeric(t(res$IND[,-1])+1),as.numeric(t(res$DIST[,-1])))
    adj<-Matrix::sparseMatrix(i=adj[1,],j=adj[2,],x=adj[3,],dims=c(N,N))
    return(adj)
}

#' Symmetrise similarity (or adjacency) matrix
#'
#' \code{knn.spadj2sym} converts a similarity matrix, as generated by, for instance, \code{knn.adj2spadjsim}, (or, optionally, an adjacency matrix) to a symmetric version thereof.
#'  For each corresponding pair of entries (row \code{i}, column \code{j} and row \code{j}, column \code{i}), we take the maximum of either entry.
#'
#' @param adj similarity or adjacency matrix (\code{sparseMatrix} object).
#'
#' @return Symmetrised matrix (\code{sparseMatrix} object).
#'
#' @examples
#' # create an artificial data set of two clusters in 2-dimensional space
#' coord <- cbind(c(rnorm(1000,100,20), rnorm(500, 20, 5)), c(rnorm(1000,100,20), rnorm(500, 120, 5)))
#'
#' sim <- knn.adj2spadjsim(     #         to similarity matrix
#'          knn.raw2adj(        #      to sparse
#'            knn.adj.raw(      #   raw adjacency matrix
#'              coord, K=5)))   # 5 nearest neighbours
#'
#' isSymmetric(sim) # FALSE
#'
#' sim <- knn.spadj2sym(sim)
#'
#' isSymmetric(sim) # TRUE
#'
#' @export
knn.spadj2sym<-function(adj){
    g1<-igraph::graph_from_adjacency_matrix(adj,weighted=TRUE,mode="max")

    adj<-igraph::as_adjacency_matrix(g1,attr="weight")
    return(adj)
}

##METHOD CHANGED
knn.spadj.symmetrize<-function(adj){
    adj<-(Matrix::t(adj)+adj)/2
}

knn.adj2spadj<-function(adj){
    N<-max(adj[1,])
    Matrix::sparseMatrix(i=adj[1,],j=adj[2,],x=adj[3,],dims=c(N,N))
}

knn.adj2maxgraph<-function(adj,mode="max"){
    N<-max(adj[1,])
    igraph::graph_from_adjacency_matrix(Matrix::sparseMatrix(i=adj[1,],j=adj[2,],x=adj[3,],dims=c(N,N)),weighted=TRUE,mode=mode)
}
#' Convert sparse adjacency matrix to similarity matrix
#'
#' \code{knn.adj2spadjsim} converts a sparse adjacency matrix, as generated by \code{knn.raw2adj}, to a similarity matrix using some kernel function.
#'  Distinctly from an adjacency matrix, entries in a similarity matrix are dependent upon, but not equivalent to, distances.
#'  They are calculated using a probability density function (PDF). Standard deviation of distribution is taken as median distance or some function thereof.
#'  The particular distance value is supplied as the quantile parameter for the PDF.
#'
#' @param adj sparse adjacency matrix, as generated by \code{knn.raw2adj}.
#' @param kernel kernel function to use for computing similarity.
#'
#' @return A ``sparseMatrix`` object containing similarities of each pair of points from the original adjacency matrix.
#'
#' @examples
#' # create an artificial data set of two clusters in 2-dimensional space
#' coord <- cbind(c(rnorm(1000,100,20), rnorm(500, 20, 5)), c(rnorm(1000,100,20), rnorm(500, 120, 5)))
#'
#' sim <- knn.spadj2sym(  #             to symmetric
#'   knn.adj2spadjsim(    #          to similarity matrix
#'    knn.raw2adj(        #       to sparse
#'      knn.adj.raw(      #    raw adjacency matrix
#'        coords, K=5))))
#'
#' @export
knn.adj2spadjsim<-function(adj,kernel=c("Exp"),epsilon=NULL){
    N<-max(adj[1,])

    if (kernel=="SE") {
        if (is.null(epsilon)) epsilon<-2*median(adj[3,])^2 ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]^2/epsilon)
    }
    ##METHOD CHANGE
    if (kernel=="SE0") {
        if (is.null(epsilon)) epsilon<-median(adj[3,]) ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]^2/epsilon)
    }
    if (kernel=="Lap") {
        if (is.null(epsilon)) epsilon<-median(adj[3,]) ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]/epsilon)
    }

    if (kernel=="Exp") {
        if (is.null(epsilon)) epsilon<-2*median(adj[3,])^2 ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]/epsilon)
    }

      if (kernel=="ExpM") {
        if (is.null(epsilon)) epsilon<-max(adj[3,]) ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]/epsilon)
      }

    if (kernel=="ExpMr") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-apply(mm,MARGIN=1,max)
        adj[3,]<-as.numeric(t(t(mm)/mmm))
        adj[3,]<-exp(-adj[3,])
    }
     ##METHOD CHANGE
    if (kernel=="ExpMer") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-apply(mm,MARGIN=1,median)
        adj[3,]<-as.numeric(t(t(mm)/mmm))
        adj[3,]<-exp(-adj[3,])
    }
    ##METHOD CHANGE
    if (kernel=="ExpMer2") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-2*apply(mm,MARGIN=1,median)
        adj[3,]<-as.numeric(t(t(mm)/mmm))
        adj[3,]<-exp(-adj[3,])
      }
     if (kernel=="SEMr") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-apply(mm,MARGIN=1,max)
        adj[3,]<-as.numeric(t(t(mm^2)/mmm))
        adj[3,]<-exp(-adj[3,])
     }

    ##METHOD CHANGE
      if (kernel=="SEMer") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-apply(mm,MARGIN=1,median)^2
        adj[3,]<-as.numeric(t(t(mm^2)/mmm))
        adj[3,]<-exp(-adj[3,])
      }

    ##METHOD CHANGE
      if (kernel=="SEMer0") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-apply(mm,MARGIN=1,median)
        adj[3,]<-as.numeric(t(t(mm)/mmm))
        adj[3,]<-exp(-adj[3,])
      }
    ##METHOD CHANGE
      if (kernel=="SEMer2") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-2*(apply(mm,MARGIN=1,median)^2)
        adj[3,]<-as.numeric(t(t(mm^2)/mmm))
        adj[3,]<-exp(-adj[3,])
      }
    if (kernel=="ExpSigma"){
         D<- matrix(adj[3,],nrow=N)
         sigmaV<-C_sigma(D,log2(ncol(D)))
         adj[3,]<-as.numeric(t(t(D)/sigmaV))
         adj[3,]<-exp(-adj[3,])
    }
    Matrix::sparseMatrix(i=adj[1,],j=adj[2,],x=adj[3,],dims=c(N,N))
}

## knn.spsim.construct<-function(XX,K,metric="euclid",epsilon=NULL){
##     D<-ncol(XX)
##     N<-nrow(XX)
##     print(N)
##     print(D)
##     res<-NULL
##     if (metric=="euclid"){ res<-.Call("_tviblindi_createKNNgraph_eu",X=as.numeric(t(XX)),D=D,K=K+1,N=N)}
##     else{
##         if (metric=="cosine") res<-.Call("_tviblindi_createKNNgraph_cos",X=as.numeric(t(XX)),D=D,K=K+1,N=N)
##     }

##     adj<-rbind(rep(1:N,each=K),as.numeric(t(res$IND[,-1])+1),as.numeric(t(res$DIST[,-1])))
##     if (is.null(epsilon)) epsilon<-2*median(adj[3,])^2 ## underestimating epsilon - sample!
##     ## epsilon=0.1
##     print(epsilon)
##     adj[3,]<-exp(-adj[3,]/epsilon)
##     adj<-sparseMatrix(i=adj[1,],j=adj[2,],x=adj[3,],dims=c(N,N))
##     return(adj)
## }

#' Generate a sparse Laplacian matrix from an adjacency matrix
#'
#' \code{sparse.Laplacian.construct} produces a generalised graph Laplacian, which can be seen as a symmetrised transition matrix. Transition probabilities are computed using a specified kernel function.
#'
#' @param adj adjacency matrix, as generated by \code{knn.raw2adj}.
#' @param kernel kernel function for calculating transition probabilities: \code{Exp} for exponential kernel, \code{SE} for SE kernel or \code{Lap} for Laplace kernel. Defaults to \code{Exp}.
#'
#' @return A \code{sparseMatrix} object representing a symmetric graph Laplacian.
#'
#' @export
sparse.Laplacian.construct<-function(adj,kernel=c("Exp"),epsilon=NULL){
    cat("constructing Laplacian\n")
    N<-max(adj[1,])

    if (kernel=="SE") {
        if (is.null(epsilon)) epsilon<-2*median(adj[3,])^2 ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]^2/epsilon)
    }

    if (kernel=="Lap") {
        if (is.null(epsilon)) epsilon<-median(adj[3,]) ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]/epsilon)
    }

     if (kernel=="Exp") {
        if (is.null(epsilon)) epsilon<-2*median(adj[3,])^2 ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]/epsilon)
    }

     if (kernel=="ExpM") {
        if (is.null(epsilon)) epsilon<-max(adj[3,]) ## underestimating epsilon - sample!
        adj[3,]<-exp(-adj[3,]/epsilon)
     }

     if (kernel=="ExpMr") {
        mm<- matrix(adj[3,],nrow=N)
        mmm<-apply(mm,MARGIN=1,max)
        adj[3,]<-as.numeric(t(t(mm)/mmm))
        adj[3,]<-exp(-adj[3,])
     }

     if (kernel=="ExpSigma"){
         D<- matrix(adj[3,],nrow=N)
         sigmaV<-C_sigma(D,log2(ncol(D)))
         adj[3,]<-as.numeric(t(t(D)/sigmaV))
         adj[3,]<-exp(-adj[3,])
    }

    cat("transition constructed\n")
    K<-Matrix::sparseMatrix(i=adj[1,],j=adj[2,],x=adj[3,],dims=c(N,N))
    rm(adj)
    cat("sparse constructed\n")
    g1<-igraph::graph_from_adjacency_matrix(K,weighted=TRUE,mode="max")
    K<-igraph::as_adj(g1,attr="weight")
    rm(g1)
    cat("symmetrized\n")
    diag(K)<-1
    v<-Matrix::rowSums(K)
    Dv<-Matrix::Diagonal(x=1/sqrt(v))
    K= Dv %*% K %*% Dv  # symmetric graph Laplacian
    ## v<-Matrix::rowSums(K)
    ## Dv<-Diagonal(x=1/v)
    ## K <- Dv %*% K #normalized laplacian
    ## would lose symmetry and this one is adjoint!
    return(K)

}



## sparse.Laplacian_alpha.construct<-function(adj,alpha=1){
##     cat("constructing Laplacian\n")
##     N<-max(adj[1,])
##     epsilon<-2*median(adj[3,])^2 ## underestimating epsilon - sample!
##     adj[3,]<-exp(-adj[3,]^2/epsilon)
##     cat("transition constructed\n")
##     K<-Matrix::sparseMatrix(i=adj[1,],j=adj[2,],x=adj[3,],dims=c(N,N))
##     rm(adj)
##     cat("sparse comstructed\n")
##     g1<-igraph::graph_from_adjacency_matrix(K,weighted=TRUE,mode="max")
##     K<-igraph::as_adj(g1,attr="weight")
##     rm(g1)
##     cat("symmetrized\n")
##     diag(K)<-1
##     v<-Matrix::rowSums(K)
##     Dv<-Matrix::Diagonal(x=1/sqrt(v))
##     K=  Matrix::Diagonal(nrow(K)) - Dv %*% K %*% Dv  # symmetric graph Laplacian
##     ## v<-Matrix::rowSums(K)
##     ## Dv<-Diagonal(x=1/v)
##     ## K <- Dv %*% K #normalized laplacian
##     ## would lose symmetry and this one is adjoint!
##     return(K)

## }

sparse.Laplacian.construct.sym<-function(K){
    v<-Matrix::rowSums(K)
    Dv<-Matrix::Diagonal(x=1/sqrt(v))
    K= Dv %*% K %*% Dv
    return(K)

}

#' Generate diffusion matrix from graph Laplacian
#'
#' \code{sparse.diffuse} computes a diffusion matrix of a Markov chain, given some generalised graph Laplacian, as generated by \code{sparse.Laplacian.construct}. Diffusion coordinates are extracted via
#' eigendecomposition of the Laplacian. Diffusion process is then simulated by raising the transition matrix to a power of \code{t} (this is to say, \code{t} 'steps'
#' in the diffusion process are simulated).
#'  Adapted from package \code{diffusionMap}.
#'
#' @param Asp generalised graph Laplacian matrix, as generated by \code{sparse.Laplacian.construct}.
#' @param neigen number of eigenvectors to be computed. Defaults to NULL, in which case the number of eigenvectors computed is determined by either \code{maxdim} + 1 (the first left eigenvector is disregarded,
#' since it represents stationary distribution of the Markov chain), or dimension of \code{Asp} (minimum of these values is taken). In simulation of the diffusion process, a value of \code{neigen} corresponding
#' to less than a 95-percent dropoff in eigenvalue versus maximum is taken.
#' @param t number of steps in the diffusion process. Defaults to 0. If t <= 0, use multi-scale geometry.
#' @param maxdim maximum dimension of diffusion coordinates.
#' @param method method employed for eigendecomposition. Currently supports \code{arpack} from package \code{igraph}, else default to \code{eigs_sym} from package \code{PRIMME}.
#'
#' @return
#' \code{sparse.diffuse} returns a list containing the following elements.
#'
#' \code{X} is a matrix of diffusion coordinates.
#'
#' \code{phi0} is the left eigenvector representing stationary distribution of the diffusion process.
#'
#' \code{eigenvals} is a vector containing eigenvalues of \code{Asp}, except for the first (largest) value, corresponding to \code{phi0}.
#'
#' \code{psi} is a matrix of right eigenvectors (by rows).
#'
#' \code{phi} is a matrix of left eigenvectors (by columns).
#'
#' \code{neigen} is the \code{neigen} input parameter.
#'
#' \code{eigenmult} are eigen-multipliers of the diffusion map.
#'
#' @export
sparse.diffuse<- function(Asp, neigen=NULL, t=0, maxdim=50,method="arpack") {
    n=dim(Asp)[1]
    start = proc.time()[3]
    ##from package diffusionMap - modified
    f = function(x, A = NULL){ # matrix multiplication for ARPACK

        as.matrix(A %*% x)
    }

    cat('Performing eigendecomposition\n') # eigendecomposition
    if(is.null(neigen)){
        neff = min(maxdim+1,n)
    }else{
        neff =  min(neigen+1, n)
    }

                                        # eigendecomposition using ARPACK
    if (method=="arpack"){
        decomp = igraph::arpack(f,extra=Asp,sym=TRUE,
                                options=list(which='LA',nev=neff,n=n,ncv=max(min(c(n,4*neff)))))
    } else {
        decomp=eigs_sym(Asp,NEig=neff)
    }

    ## lastone<-max(which(decomp$values==1.0))

    ## if (lastone>1) {
    ##     decomp$vectors<-decomp$vector[,-c(2:lastone)]
    ##     decomp$values<-decomp$values[,-c(2:lastone)]
    ##     neigen<-length(decomp$values)-1
    ## }

    psi = decomp$vectors/(decomp$vectors[,1]%*%matrix(1,1,neff))#right ev
    phi = decomp$vectors * (decomp$vectors[,1]%*%matrix(1,1,neff))#left ev
    eigenvals = decomp$values #eigenvalues

    cat('Computing Diffusion Coordinates\n')
    if(t<=0){# use multi-scale geometry
        lambda=eigenvals[-1]/(1-eigenvals[-1])
        lambda=rep(1,n)%*%t(lambda)
        if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
            lam = lambda[1,]/lambda[1,1]
            neigen = min(which(lam<.05)) # default number of eigenvalues
            neigen = min(neigen,maxdim)
            eigenvals = eigenvals[1:(neigen+1)]
            cat('Used default value:',neigen,'dimensions\n')
        }
        X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
    }
    else{# use fixed scale t
        lambda=eigenvals[-1]^t
        lambda=rep(1,n)%*%t(lambda)

        if(is.null(neigen)){#use no. of dimensions corresponding to 95% dropoff
            lam = lambda[1,]/lambda[1,1]
            neigen = min(which(lam<.05)) # default number of eigenvalues
            neigen = min(neigen,maxdim)
            eigenvals = eigenvals[1:(neigen+1)]
            cat('Used default value:',neigen,'dimensions\n')
        }

        X = psi[,2:(neigen+1)]*lambda[,1:neigen] #diffusion coords. X
    }
    cat('Elapsed time:',signif(proc.time()[3]-start,digits=4),'seconds\n')

    y = list(X=X,phi0=phi[,1],eigenvals=eigenvals[-1],eigenmult=lambda[1,1:neigen],
             psi=psi,phi=phi,neigen=neigen)
    class(y) = "dmap"
    return(y)


}

#' Generate diffusion matrix from coordinates
#'
#' \code{knn.diffuseX} generates a diffusion matrix. Diffusion coordinates are extracted via eigendecomposition of a generalised Laplacian of a K-nearest-neighbours graph.
#'  Adapted from package \code{diffusionMap}.
#'
#' @param X matrix of coordinates. Rows correspond to vertices.
#' @param K number of neighbours to use in construction a kNN graph.
#' @param neigen number of eigenvectors and eigenvalues to be computed. Defaults to number of dimensions.
#' @param metric distance metric. \code{euclid} or \code{cosine}.
#'
#' @return
#' \code{sparse.diffuse} returns a list containing the following elements.
#'
#' \code{X} is a matrix of diffusion coordinates.
#'
#' \code{phi0} is the left eigenvector representing
#' stationary distribution of the diffusion process.
#'
#' \code{eigenvals} is a vector containing eigenvalues of \code{Asp}, except for the first (largest) value, corresponding to \code{phi0}.
#'
#' \code{psi} is a matrix of right eigenvectors as rows.
#'
#' \code{phi} is a matrix of left eigenvectors as columns.
#'
#' \code{neigen} is the input parameter \code{neigen}.
#'
#' \code{eigenmult} are eigen-multipliers of the diffusion map.
#'
#' @export
knn.diffuseX<-function(X,K=30,neigen=dim(X)[2],metric="euclid",t=0){
    sparse.diffuse(sparse.Laplacian.construct(knn.raw2adj(knn.adj.raw.parallel(X,K,metric=metric))),neigen=neigen,t=t)
}

#' Convert sparse similarity matrix to "triples" matrix of edges
#'
#' \code{sparse2triples} converts a \code{sparseMatrix} object (typically, a similarity matrix) to a "triples" format: a matrix describing oriented graph edges.
#'  Columns correspond to vertices and edge weight.
#'
#' @param m similarity (or analogous) matrix (\code{sparseMatrix} object).
#'
#' @return \code{sparse2triples} generates a 3-column matrix. Columns correspond to \code{from} vertex, \code{to} vertex and edge \code{weight}, respectively.
#'
#' @examples
#' # create an artificial data set of two clusters in 2-dimensional space
#' coord <- make_snowman(1000)
#' origin <- 1
#'
#' # create a similarity matrix
#' sim <- knn.spadj2sym(  #             to symmetric
#'   knn.adj2spadjsim(    #          to similarity matrix
#'    knn.raw2adj(        #       to sparse
#'      knn.adj.raw(      #    raw adjacency matrix
#'        coords, K=5))))
#'
#' pseudotime <- assign_distance(sim, origin) # assign pseudotime
#' oriented <- remove_back_R(sparse2triples(sim), # change matrix format
#'                             pseudotime$res)      # orient the graph according to pseudotime values
#'
#' @export
sparse2triples <- function(m) {
    m<-Matrix::summary(m)
    colnames(m)<-c("from","to","weight")
    m
}

#' Remove backward edges from "triples" matrix
#'
#' \code{remove_back_R} takes a "triples" matrix of edges, as generated by \code{sparse2triples} and \code{time}, as generated by \code{assign_distance}.
#'  Edges such that the \code{from} vertex is assigned a distance value by \code{time} lower or equal to value assigned to the \code{to} vertex are removed.
#'
#' @param triplets edge matrix, as generated by \code{sparse2triples}.
#' @param time distance assignment, as generated by \code{assign_distance}.
#'
#' @return \code{remove_back_R} outputs the \code{triplets} matrix, possibly with some rows removed.
#'
#' @examples
#' # create an artificial data set of two clusters in 2-dimensional space
#' coord <- make_snowman(1000)
#' origin <- 1
#'
#' # create a similarity matrix
#' sim <- knn.spadj2sym(  #             to symmetric
#'   knn.adj2spadjsim(    #          to similarity matrix
#'    knn.raw2adj(        #       to sparse
#'      knn.adj.raw(      #    raw adjacency matrix
#'        coords, K=5))))
#'
#' pseudotime <- assign_distance(sim, origin) # assign pseudotime
#' oriented <- remove_back_R(sparse2triples(sim), # change matrix format
#'                             pseudotime$res)      # orient the graph according to pseudotime values
#'
#' @export
remove_back_R<-function(triplets,time){
    triplets[which(time[triplets[,1]]<time[triplets[,2]]),]
}

KofRawN<-function(rawadj,K=NULL){
    if (is.null(K)) K<-ncol(rawadj$IND)
    if ((K+1) >= ncol(rawadj$IND)){
        warning("K+1 >= N")
        K<-ncol(rawadj$IND)-1
    }
    rawadj$IND<-rawadj$IND[,1:(K+1)]
    rawadj$DIST<-rawadj$DIST[,1:(K+1)]
    return(rawadj)
}
