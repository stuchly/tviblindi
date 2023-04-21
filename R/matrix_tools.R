orient.sim.matrix<-function(sim,pseudotime,breaks=NULL,base=2){
  N<-length(pseudotime$res)
  oriented <- remove_back_R(sparse2triples(sim), pseudotime$res)
  if (!is.null(breaks)){
    breaks<-as.numeric(cut(as.numeric(as.factor(pseudotime$res)),breaks))
    breaks<-breaks[oriented$to]-breaks[oriented$from]
    oriented$weight<-oriented$weight/base^breaks
  }
  A<-Matrix::sparseMatrix(i=oriented$from, j=oriented$to, x=oriented$weight, dims=c(N,N))
  #is it ok?
  .DD<-Matrix::rowSums(A)
  .DD <- Matrix::Diagonal(x=.DD^-1)
  A<-.DD%*%A
}

transition.matrix<-function (x, K = 30, kernel = "SEMer", sym = "max",kepsilon = NULL,normalize=TRUE)
{

  if (K > dim(x$KNN$IND)[2]) {
    K <- min(K, dim(x$KNN)[2])
    warning("K > dim(KNN)[2]; K<-min(K,dim(x$KNN)[2])")
  }
  d <- KofRawN(x$KNN, K)
  d <- knn.raw2adj(d)
  dsym <- knn.spadj2sym(knn.adj2spadj(d))
  symB <- TRUE
  if (sym == "none") {
    sim <- knn.adj2spadjsim(d, kernel = kernel, epsilon = kepsilon)
    symB <- FALSE
  }
  else if (sym == "mean")
    sim <- knn.spadj.symmetrize(knn.adj2spadjsim(d, kernel = kernel,
                                                 epsilon = kepsilon))
  else if (sym == "prob")
    sim <- knn.spadj.symmetrize.P(knn.adj2spadjsim(d, kernel = kernel,
                                                   epsilon = kepsilon))
  else if (sym == "max")
    sim <- knn.spadj2sym(knn.adj2spadjsim(d, kernel = kernel,
                                          epsilon = kepsilon))
  else if (sym == "min") {
    d <- t(summary(dsym))
    sim <- knn.adj2spadjsim1(d, kernel = kernel, epsilon = kepsilon)
    symB = FALSE
  }
  else stop("symmetrisation not implemented")
  if (normalize){
    diag(sim)<-1
    .DD <- Matrix::rowSums(sim)

    .DD <- Matrix::Diagonal(x=.DD^-1)
    return(.DD%*%sim)
  } else return(sim)
}
