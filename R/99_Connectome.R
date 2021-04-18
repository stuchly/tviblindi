#' Create connectome-analysis object from a tviblindi object
#'
#' Generates a diagram of connections between Louvain clusters of high-dimensional expression data, indicating flow (and direction) between tentative developmental stages.
#' This is based on previously simulated random walks (stochastic traversals of a k-NN graph governed by properties of a previously generated transition model).
#' If \code{equinumerous} was set to true in random walks (equal number of walks per terminus enforced), this biases the results (undesirable).
#'
#' @param tv trajectory-inference object of class \code{tviblindi}
#' @param k integer: number of nearest neighbors for Louvain clustering
#' @param layout_name string: name of 2-dimensional layout to use for connectome analysis
#' @param labels_name string: name of labels vector used to display pie charts for each connectome graph node
#' @param transition_model_name string: name of transition model to use for connectome analysis
#' @param png optional string: file path to a PNG file of the diagram to be generated. Default value is \code{NULL}, whereby the currently active plotting device is used instead
#'
#' @export
Connectome <- function(
  tv,
  k = 30,
  layout_name = 'default',
  labels_name = 'default',
  transition_model_name = 'default',
  png = NULL
) {
  if (is.null(tv$metaclusters)) { 
    .msg('Louvain clustering...', endline = FALSE)
    if (is.null(tv$dsym)) {
      d <- KofRawN(tv$kNN, k)
      d <- rbind(rep(seq_len(nrow(d$idcs)), each = k), as.numeric(t(d$idcs[, -1]) + 1), as.numeric(t(d$dist[, -1])))
      d <- Matrix::sparseMatrix(i = d[1, ], j = d[2, ], x = d[3, ], dims = rep(max(d[1, ]), times = 2))
      dsym <- igraph::as_adjacency_matrix(
        igraph::graph_from_adjacency_matrix(d, weighted = TRUE, mode = 'max'),
        attr = 'weight'
      )
    } else {
      tv$dsym <- dsym
    }
    g <- igraph::graph_from_adjacency_matrix(dsym, weighted=TRUE, mode = 'undirected')
    tv$metaclusters <- igraph::cluster_louvain(g)$membership
    .msg_alt_good(' done!')
    rm(g)
  } else {
    dsym <- tv$dsym
  }
  
  ends <- c(tv$walks[[transition_model_name]]$starts[-1] - 1, length(tv$walks[[transition_model_name]]$v))
  starts <- tv$walks[[transition_model_name]]$starts
  
  edge_list <- 
    rbind(
      cbind(
        tv$walks[[transition_model_name]]$v[-ends],
        tv$walks[[transition_model_name]]$v[-starts]
      ),
      dim(dsym)
    )
  g <- igraph::graph_from_edgelist(edge_list, directed = TRUE)
  g <- igraph::contract(g, tv$metaclusters)
  g <- igraph::as_adjacency_matrix(g)
  s <- Matrix::summary(g)
  s <- s[s$i != s$j, ]
  
  clus <- sort(unique(tv$metaclusters))
  
  s <- Matrix::sparseMatrix(i = s$i, j = s$j, x = s$x, dims = rep(length(clus), times = 2))
  
  g_layout <- do.call(rbind, lapply(clus, function(idx) colMeans(tv$layout[[layout_name]][which(tv$metaclusters == idx), ])))
  
  a <- Matrix::Diagonal(x = Matrix::rowSums(Matrix::t(s) + s)^(-1)) %*% s
  gg <- igraph::graph_from_adjacency_matrix(a, weighted = TRUE, mode = 'directed')
  igraph::E(gg)$width <- igraph::E(gg)$weight^1.2 * 19
  igraph::E(gg)$curved <- TRUE
  igraph::E(gg)$arrow.size <- igraph::E(gg)$width + 0.2
  igraph::E(gg)$arrow.width <- 0.7
  igraph::V(gg)$label.cex <- 1
  igraph::V(gg)$label.dist <- 0.8
  igraph::V(gg)$label.color <- 'black'
  igraph::V(gg)$label.family <- 'sans'
  
  pies <- lapply(clus, function(idx) as.vector(table(tv$labels[[labels_name]][tv$metaclusters == idx])))
  
  colours <- rep(list(rainbow(nlevels(tv$labels[[labels_name]]) + 1)), length(clus))
  
  lsize <- 0.4
  if (!is.null(png)) {
    png(png, 2000, 2000)
    lsize <- 1.5
  }
  old_mar <- par(no.readonly = TRUE)$mar
  par(mar = c(0, 0, 2, 0))
  igraph::plot.igraph(
    gg, layout = g_layout, main = 'Inferred conectome', vertex.size = 7, vertex.shape = 'pie', vertex.pie = pies,
    vertex.pie.color = colours
  )
  par(mar = old_mar)
  legend('topright', legend = levels(tv$labels[[labels_name]]), col = colours[[1]], pch = 19, cex = lsize)
  if (!is.null(png))
    dev.off()
  
  invisible(tv)
}

KofRawN<-function(
  rawadj, K = NULL
) {
  if (is.null(K))
    K <- ncol(rawadj$idcs)
  if ((K+1) >= ncol(rawadj$idcs)){
    warning("K+1 >= N")
    K <- ncol(rawadj$idcs) - 1
  }
  rawadj$idcs <- rawadj$idcs[, 1:(K + 1)]
  rawadj$dist <- rawadj$dist[, 1:(K + 1)]
  rawadj
}
