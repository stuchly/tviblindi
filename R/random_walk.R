random_walk_adj_R<-function(adj,start,steps){
    vv<-c(start,rep(NA,steps))
    step<-1
    while(step<=steps){
        step<-step+1
        x<-adj[start,]
        ss<-which(x!=0)
        if (length(ss)==0) break
        start<-ss[sample(x=1:length(ss),size=1,prob = as.numeric(x[ss]))]
        vv[step]<-start
    }
    return(vv[!is.na(vv)])
}

random_walk_adj<-function(adj,start,steps){
    C_random_walk_adj(Matrix::t(adj),start,steps)$v
}

random_walk_adj_u<-function(adj,start,steps){
    C_random_walk_adj_u(Matrix::t(adj),start,steps)$v
}

random_walk_adj_N<-function(adj,start,steps,N){
    C_random_walk_adj_N(Matrix::t(adj),start,steps,N)
}

random_walk_adj_N_push<-function(adj,start,N){
    C_random_walk_adj_N_push(Matrix::t(adj),start,N)
}

random_walk_adj_N_push_std<-function(adj,start,N){
    C_random_walk_adj_N_push_std(Matrix::t(adj),start,N)
}

#' Generate random walks through graph
#'
#' \code{random_walk_N} generates a set of random walks through a graph, provided an oriented transition matrix, index of origin vertex (starting point
#'  for all paths) and number of walks to simulate. Since we are working with a connected oriented graph, there exist end-vertices where a walk must
#'  terminate. Therefore, we simply use entries in the transition matrix to calculate probability of transitioning between vertices and simulate a pre-set
#'  number of walks iteratively.
#'
#' @param A transition matrix.
#' @param origin index of the origin vertex.
#' @param N number of walks to simulate.
#'
#' @return \code{random_walk_N} returns a list containing the following elements. \code{v} is an integer vector of all simulated walks in terms of vertex
#'  indices. \code{starts} is an integer vector containing indices of \code{v} which correspond to start-points of each simulated walk.
#'
#' @export
random_walk_N<-function(A,origin,N){
  C_random_walk(Matrix::t(A),origin,N)
}
