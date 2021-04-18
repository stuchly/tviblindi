.update_walks_by_termini <- function(
  tv,
  pathmodel_name,
  pseudotime,
  marked_termini,
  termini_per_path,
  death_birth_ratio,
  death_on_x_axis
) {
  ## Identify chosen walks
  idcs <- which(termini_per_path %in% marked_termini)
  
  ## Remove other walks
  walks.selected <- lapply(idcs, function(idx) select_walks_nodes(tv$walks[[pathmodel_name]], idx))
  lens <- sapply(walks.selected, length)
  walks.selected <- list(v = unlist(walks.selected), starts = c(1, 1 + cumsum(lens[-length(lens)])))
  
  ## Put terminal node with highest pseudotime at end of each walk
  if (length(marked_termini) > 1) {
    max_t <- marked_termini[which.max(pseudotime$res[marked_termini])][1]
    ends <- c(walks.selected$starts[-1] - 1, length(walks.selected$v))
    walks.selected$v[ends] <- max_t
  }
  
  ## Cluster and triangulate walks
  withProgress(message = 'Contracting trajectories', expr = {
    walks_clusters <- remove_cycles(contract_walks(walks.selected, tv$clusters))
    
    sel <- 1:length(walks_clusters$starts)
    s2 <- which(unlist(lapply(tv$filtration$cmplx, function(x) length(x) == 2))) # 2-simplex idcs
    cmplx <- tv$filtration$cmplx[s2] # 2-simplices
    cmplx_hashed <- hash_cmplx(cmplx) # hash the simplices list for faster look-up
    edges <- do.call(rbind, lapply(cmplx, function(sim, XX) c(sim, stats::dist(tv$codes[sim,])), XX = tv$codes))
    colnames(edges) <-c('from', 'to', 'weight')
    
    graph <- igraph::graph_from_edgelist(edges[, 1:2])
    igraph::E(graph)$weight <- edges[, 3]
    triangulation <- list()[1:length(sel)]
    
    j <- 0
  })
  
  withProgress(message = 'Triangulating paths', expr = {
    for (i in sel) {
      j <- j + 1
      triangulation[[j]] <-
        s2[triangulate_trajectory(select_walks_nodes(walks_clusters, i), tv$codes, cmplx_hashed, graph)]
    }
  })
  
  N <- length(idcs)
  repre <- list()[1:N]
  repre[[1]] <- integer(0)
  j <- 1
  
  N <- length(idcs)
  tick <- 1 / N
  
  withProgress(message = 'Computing representations', value = tick, expr = {
    for (idx in 2:N) {
      j <- j + 1
      cycle <- path_sum(triangulation[[1]], triangulation[[idx]])
      this_repre <- get_rep_straight(cycle, tv$reduced_boundary, tv$boundary, update = TRUE)
      repre[[j]] <- this_repre
      incProgress(tick)
    }
  })
  
  walks <- lapply(1:N, function(idx) { select_walks_nodes(walks.selected, idx) })
  
  p <- .compute_persistence(tv, repre, death_birth_ratio, death_on_x_axis)
  
  list(
    random_walks  = walks,
    repre = repre,
    pers = p$pers,
    pers_diag = p$pd,
    triangulation = triangulation
  )
}

.update_walks_by_dendrogram <-  function(
  tv,
  pathmodel_name,
  pseudotime,
  marked_termini,
  idcs,
  walks_sel,
  triangulation,
  termini_per_path,
  death_birth_ratio,
  death_on_x_axis
) {
  ## Identify chosen walks
  triangulation<-triangulation[idcs]
  
  idcs <- walks_sel[idcs]
  
  ## Remove other walks
  walks.selected <- lapply(idcs, function(idx) select_walks_nodes(tv$walks[[pathmodel_name]], idx))
  lens <- sapply(walks.selected, length)
  walks.selected <- list(v      = unlist(walks.selected),
                         starts = c(1, 1 + cumsum(lens[-length(lens)])))
  
  ## Put terminal node with highest pseudotime at end of each walk
  if (length(marked_termini) > 1) {
    max_t <- marked_termini[which.max(pseudotime$res[marked_termini])][1]
    ends <- c(walks.selected$starts[-1] - 1, length(walks.selected$v))
    walks.selected$v[ends] <- max_t
  }
  
  N <- length(idcs)
  repre <- list()[1:N]
  repre[[1]] <- integer(0)
  j <- 1
  
  N <- length(idcs)
  tick <- 1 / N
  
  withProgress(message = 'Computing representations', value = tick, expr = {
    for (idx in 2:N) {
      j          <- j + 1
      cycle      <- path_sum(triangulation[[1]], triangulation[[idx]])
      this_repre <- get_rep_straight(cycle, tv$reduced_boundary, tv$boundary, update = TRUE)
      repre[[j]] <- this_repre
      incProgress(tick)
    }
  })
  
  N <- length(idcs)
  walks <- lapply(1:N, function(idx) { select_walks_nodes(walks.selected, idx) })
  
  p <- .compute_persistence(tv, repre, death_birth_ratio, death_on_x_axis)
  
  list(
    random_walks = walks,
    repre = repre,
    pers = p$pers,
    pers_diag = p$pd,
    walks_sel = idcs,
    triangulation = triangulation
  )
}


.persistence_diagram <- function(
  pers,
  death_birth_ratio,
  death_on_x_axis
) {
  pd <- data.frame(
    Dimension = pers$vals$dim,
    Birth = pers$vals$birth,
    Death = pers$vals$death,
    BirthSimplex = pers$idcs$birth,
    DeathSimplex = pers$idcs$death,
    xplot = if (death_on_x_axis) { pers$vals$death } else { (pers$vals$birth + pers$vals$death) / 2 },
    yplot = if (death_birth_ratio) { pers$vals$death / pers$vals$birth } else { (pers$vals$death - pers$vals$birth) / 2 }
  )
  pd[pd$Dimension > 0, ]
}

.compute_persistence <- function(
  tv,
  repre,
  death_birth_ratio,
  death_on_x_axis
) {
  repre <-
    if (!is.null(repre))
      unique(unlist(repre))
    else
      tv$reduced_boundary$nonzero_col
  
  s <- which(
    (tv$reduced_boundary$values[tv$reduced_boundary$low] != tv$reduced_boundary$values[tv$reduced_boundary$nonzero_col]) &
     (tv$reduced_boundary$nonzero_col %in% repre)
  )
  pers <- list(
    idcs = data.frame(
      dim = tv$reduced_boundary$dim[s],
      birth = tv$reduced_boundary$low[s],
      death = tv$reduced_boundary$nonzero_col[s]
    ),
    vals = data.frame(
      dim = tv$reduced_boundary$dim[s],
      birth = tv$reduced_boundary$values[tv$reduced_boundary$low[s]],
      death = tv$reduced_boundary$values[tv$reduced_boundary$nonzero_col[s]]
    )
  )
  pd <- .persistence_diagram(pers, death_birth_ratio, death_on_x_axis)
  list(
    pers = pers,
    pd = pd
  )
}


#' Generate triangulated trajectory
#'
#' \code{triangulate_trajectory} re-creates a random walk usign 1-simplices generated by a filtration function, Dijkstra's shortest paths search.
#' The function outputs a 1-chain which approximates and represent their topology.
#' Since all vertices from the original walk must be included, the triangulated trajectory must be of the same or of greater length.
#'
#' @param points integer vector: indices of the nodes which comprise a random walk (vertices of the original point cloud)
#' @param cluster_codes numeric matrix: cluster-center coordinates, with each row corresponding to a single cluster center
#' @param cmplx_hashed a hashed environment with the simplices in a filtration, as generated by \code{tviblindi:::hash_cmplx}.
#' @param graph an \code{igraph} object representing the graph in which to run the shortest paths search (implemented in \code{tviblindi:::.update_walks_by_termini})
#'
#' @return integer vector of graph node indices
#'
#' @references
#' \insertRef{Csardi2006}{tviblindi}
#' 
#' \insertRef{Eddelbuettel2013}{tviblindi}
#'
#' @export
triangulate_trajectory <- function(
  points,
  cluster_codes,
  cmplx_hashed,
  graph
) {
  keys <- paste(points[-length(points)], points[-1], sep = '_')
  values <- cbind(points[-length(points)], points[-1])
  res <- NULL

  for (i in 1:nrow(values)) {
    idx <- cmplx_hashed[[keys[i]]]
    if (!is.null(idx)) {
      res <- c(res, idx)
      next
    }

    idx <- as.integer(
      igraph::shortest_paths(
        graph, from = values[i, 1], to = values[i, 2], mode = 'all', output = 'epath'
      )$epath[[1]]
    )
    res <- c(res, idx)
  }

  res.t <- table(res)

  cycles <- which(res.t %% 2 == 0)
  if (length(cycles) > 0) {
    cycles <- which(res %in% as.integer(names(res.t[cycles])))
    res <- res[-cycles]
  }

  sort(unique(res))
}


#' Contract random walks using result of clustering
#'
#' \code{contract_walks} re-creates random walks through a graph using assignments of nodes to clusters.
#'
#' @param walks list of random walks, consisting of node indices in slot \code{v} and start indices referring to entries from \code{v} in \code{starts}
#' @param clusters integer vector: assignments of each node along the walk to a cluster
#'
#' @return modified \code{walks} list
#'
#' @export
contract_walks <- function(
  walks,
  clusters
) {
  contracted <- list()
  
  for (i in 1:length(walks$starts)) {
    w     <- select_walks_nodes(walks, i)
    edges <- lapply(1:(length(w) - 1), function(j, k = j + 1) clusters[w[j:k]]) # walk as list of edges (using cluster numbering)
    edges <- edges[unlist(lapply(edges, function(j) j[1] != j[2]))]             # remove edges leading to self
    if (length(edges) == 0) {
      stop('Simulated walks are confined to one SOM cluster. Either the SOM grid dimension is very small or data is extremely homogeneous')
    }
    w <- c(
      edges[[1]][1],
      unlist(lapply(edges, function(i) i[2]))
    )
    contracted[[i]] <- w
  }
  
  starts <- c(1, unlist(lapply(contracted, length)))
  starts <- starts[-length(starts)]
  starts <- sapply(1:length(starts), function(i) sum(starts[1:i]))
  
  list(v = unlist(contracted), starts = starts)
}

#' Select nodes of random walks by walk indices
#'
#' \code{select_walks_nodes} returns a modified \code{walks} list, containing only those walks which are selected by their indices in the original \code{walks} object.
#'
#' @param walks list of random walks, consisting of node indices in slot \code{v} and start indices referring to entries from \code{v} in \code{starts}
#' @param selection integer vector: walk indices, with entries from \code{1:length(walks$starts)}.
#'
#' @return modified \code{walks} list
#'
#' @export
select_walks_nodes <- function(
  walks,
  selection
) {
  N <- length(walks$starts)
  start_node_idx <- function(walks, idx)  walks$starts[idx]
  end_node_idx <- function(walks, idx) if (idx == N) { length(walks$v) } else { walks$starts[idx + 1] - 1 }
  nodes_idcs <- function(walks, idx) start_node_idx(walks, idx):end_node_idx(walks, idx)
  
  unlist(
    lapply(
      selection,
      function(walk_idx)
        walks$v[nodes_idcs(walks, walk_idx)]
    )
  )
}

hash_cmplx <- function(
  cmplx
) {
  hash <- new.env(hash = TRUE)
  for (i in 1:length(cmplx)) {
    
    if (length(cmplx[[i]]) == 2) {
      assign(paste(cmplx[[i]], collapse = '_'), i, envir = hash)
      assign(paste(rev(cmplx[[i]]), collapse = '_'), i, envir = hash)
    }
  }
  hash$Nsimplices <- length(cmplx)
  hash
}