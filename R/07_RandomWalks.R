#' Simulate random walks
#'
#' This function adds new slot \code{walks} to an object of class \code{tviblindi}.
#' Random walks are simulated using a pseudotime values, generated using a previously generated transition model.
#'
#' @param tv \code{tviblindi}-class object
#' @param transition_model_name string: name of transition model to use for generating random walks. Default value is '\code{default}'
#' @param n_walks integer: number of random walks to simulate. Default value is \code{1000}
#' @param equinumerous logical: indicates whether to simulate an equal number of walks for each potential terminal node (\code{N} for each). Default value is \code{FALSE}, whereby the frequency with which walks end in any given terminal node is probabilistic (depending on the transition model)
#' @param enforced_terminal_nodes optional string vector: names of labelled populations for which the node closest to centroid of that population is forced to become a terminal node. If this is specified, these enforced terminal nodes will be the only terminal nodes. Default value is \code{NULL}
#' @param labels_name string: name of labels vector to use if parameter \code{enforced_terminal_nodes} is specified. Default value is '\code{default}'
#' @param add logical: indicates whether to add the simulated walks to those that have been generated before. Otherwise, any existing walks for the given transition model are replaced. Default value is \code{FALSE}
#'
#' @details
#' This method simulates random walks on a directed graph (i.e. only edges along which pseudotime values increase are considered).
#' If the parameter \code{equinumerous} is set to \code{TRUE}, the potential terminal nodes of walks (i.e. those edges which do not lead to any other edge in the directed graph) are identified and the following procedure is followed: 
#' for each such terminal node, the sub-component of the graph from which this end could be reached is used for simulation separately, with the same number of walks simulated for each such node.
#' The same procedure is followed if \code{enforced_terminal_nodes} is non-\code{NULL}.
#'
#' @export
SimulateRandomWalks.tviblindi <- function(
  tv,
  transition_model_name = 'default',
  n_walks = 1000,
  equinumerous = FALSE,
  enforced_terminal_nodes = NULL,
  labels_name = 'default',
  add = FALSE
) {
  stopifnot(length(tv$origin) > 0)
  stopifnot(!is.null(tv$trans))
  
  add_walks <- function(tv, walks, name) {
    if (is.null(tv$walks)) {
      tv$walks <- list()
      tv$walks[[name]] <- list(starts = NULL, v = NULL)
    }
    tv$walks[[name]]$starts <- c(tv$walks[[name]]$starts, walks$starts + length(tv$walks[[name]]$v))
    tv$walks[[name]]$v <- c(tv$walks[[name]]$v, walks$v)
    return(invisible(tv))
  }
  
  if (!equinumerous & is.null(enforced_terminal_nodes)) {
    walks <- random_walks(tv$trans[[transition_model_name]], tv$origin[[tv$trans_name_origin[[transition_model_name]]]], n_walks)
    if (add) { 
      add_walks(tv, walks, transition_model_name)
    } else {
      if (is.null(tv$walks)) tv$walks <- list()
      tv$walks[[transition_model_name]] <- walks
    }
  }
  
  fates <- NULL
  if (!is.null(enforced_terminal_nodes)) {
    if (is.character(enforced_terminal_nodes)) {
      enforced_terminal_nodes <- enforced_terminal_nodes[enforced_terminal_nodes %in% levels(tv$labels[[labels_name]])]
      if (length(enforced_terminal_nodes) == 0) {
        stop('Populations specified in enforced_terminal_nodes not found in data')
      }
      ends <- NULL
      for (end in enforced_terminal_nodes) {
        idcs <- which(tv$labels == end)
        ends <- c(ends, idcs[which.min(rowSums(t(t(tv$data[idcs, ]) - colMeans(tv$data[idcs, ]))^2))])
      }
      enforced_terminal_nodes <- ends
    }
    fates        <- enforced_terminal_nodes
    equinumerous <- TRUE
  }
  
  if (is.null(tv$fates))
    tv$fates <- list()
  
  tv$fates[[transition_model_name]] <- 
    if (is.null(fates))
      which(!(1:nrow(tv$data) %in% Matrix::summary(tv$trans[[transition_model_name]])$i)) # which nodes 'do not lead to anywhere else'
    else
      fates
  
  if (equinumerous) {
    g <- igraph::graph_from_adjacency_matrix(tv$trans[[transition_model_name]], weighted = TRUE, mode = 'directed')
    igraph::V(g)$names <- 1:nrow(tv$data)
    
    if (!add)
      tv$walks[[transition_model_name]] <- NULL
    count <- 1
    n_fates <- length(tv$fates)
    for (fate in tv$fates) {
      .msg(paste0('Fate ', count, ' out of ', n_fates))
      count <- count + 1
      gf <- igraph::subcomponent(g, fate, 'in')
      gf <- igraph::induced_subgraph(g, gf)
      nnf <- igraph::V(gf)$names
      gf <- Matrix::summary(igraph::as_adjacency_matrix(gf, attr = 'weight'))
      gf[, 1] <- nnf[gf[, 1]]
      gf[, 2] <- nnf[gf[, 2]]
      adjacencies <- Matrix::sparseMatrix(
        i = gf$i,
        j = gf$j,
        x = gf$x,
        dims = dim(tv$trans)
      )
      walks <- random_walks(adjacencies, tv$origin, n_walks)
      add_walks(tv, walks, transition_model_name)
    }
  }
  
  .msg(paste0('Simulated ', if (!equinumerous) { n_walks } else { n_walks * n_fates }, ' random walks'))
  
  if (is.null(tv$n_walks))
    tv$n_walks <- list()
  if (is.null(tv$walks_equinumerous))
    tv$walks_equinumerous <- list()
  tv$n_walks[[transition_model_name]] <- n_walks
  tv$walks_equinumerous[[transition_model_name]] <- equinumerous
  
  gc(verbose = FALSE)
  
  invisible(tv)
}

SimulateRandomWalks <- function(tv, transition_model_name, n_walks, equinumerous, enforced_terminal_nodes, labels_name, add) UseMethod('SimulateRandomWalks', tv)

#' Simulate random walks
#'
#' Given a set of oriented edges associated with transition probabilities, \code{random_walks} simulates random walks through a graph as 1st-order Markov processes.
#' You do not need to call this method directly.
#'
#' @param transition_matrix a SparseMatrix object; matrix of transition probabilities between indexed nodes of a graph.
#' @param origin index of the node which is taken as origin of each walk.
#' @param N number of walks to simulate.
#'
#' @return a list with the following slots.
#' \code{v} is an integer vector with the sequences of node indices comprising each walk, concatenated back-to-back.
#' \code{starts} is an integer vector of indices in \code{v} that denotes which entries in \code{v} correspond to the first nodes in walks.
#'
#' @export
random_walks <- function(
  transition_matrix,
  origin,
  N
) {
  C_random_walk_adj_N(Matrix::t(transition_matrix), origin, 1000, N)
}

#' Remove cycles from simulated random walks
#'
#' \code{remove_cycles} requires a list with the folowing slots: \code{v} is an integer vector with the sequences of node indices comprising each walk,
#' concatenated back-to-back; \code{starts} is an integer vector of indices in \code{v} that denotes which entries in \code{v} correspond to the first nodes in walks.
#' The function removes any cycles within the walks.
#' You do not need to call this function directly.
#'
#' @param walks list of random walks, consisting of node indices in slot \code{v} and start indices referring to entries from \code{v} in \code{starts}
#'
#' @return modified \code{walks} without cycles
#' 
#' @references
#' \insertRef{Eddelbuettel2014}{tviblindi}
#' 
#' \insertRef{Bates2014}{tviblindi}
#' 
#' \insertRef{Eddelbuettel2013}{tviblindi}
#'
#' @export
remove_cycles <- function(
  walks
) {
  ends <- c(walks$starts[-1]-1, length(walks$v))
  walks.list <- lapply(1:length(walks$starts), function(i, v = walks$v, starts = walks$starts) {
    v[starts[i]:ends[i]]
  })
  
  w <- remove_cycles_int_list(walks.list, lapply(walks.list, unique))
  
  starts <- c(1, unlist(lapply(w, length)))
  starts <- starts[-length(starts)]
  starts <- sapply(1:length(starts), function(i) sum(starts[1:i]))
  
  list(
    v = unlist(w),
    starts = starts
  )
}

path_sum <- function(
  path1,
  path2
) symdiff(path1, path2)

symdiff  <- function(
  x,
  y
) setdiff(union(x, y), intersect(x, y))
