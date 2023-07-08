#' Create connectome from tvilblindi class object
#'
#' \code{Connectome}
#' @param x tviblindi class object.
#' @param png chracter; path to the png file to plot connectome. If NULL, connectome is plotted directly. Optimised for png
#' @param K integer (default K=30); number of nearest neighbors for louvain clustering
#' @param origin_name integer/character (default 1); path model to use
#' @param labels integer/character (default 1); labels to use
#' @param layout integer/character (default 1); layout to use
#' @param clusters integer vector (default NULL); provide clustering externally
#' @param arrow.sizefactor numeric vector (default 1); adjust size of an arrow
#' @param legend.position (default "topright")
#' @param lcex numeric( default 1); legend cex
#' @param notplot character (character vector) (default "ungated"); which populations should not be shown in piecharts
#' @param qq numeric (default 0); only arrows of edges width above qq percentile will be plotted
#' @param directed boolean (default TRUE); plot only edges of major direction
#' @param cluster_labels boolean; print the cluster labels (numbers) into piecharts
#' @details Computes louvain clusters and estimates the flow (and the direction) between them from simulated walks
#' (parameter \code{equinumerous} in \code{Walks}) would bias the result!
#'
#' @return \code{tviblindi} returns an invisible tviblindi class object.
#'
#' @export
connectome <-
  function (x,
            png = "connectome.png",
            K = 30,
            origin_name = 1,
            labels = 1,
            layout = 1,
            clusters = NULL,
            arrow.sizefactor = 1,
            legend.position = "topright",
            lcex = 1,
            notplot = "ungated",
            qq = 0,
            directed = TRUE,
            scale = "max",
            hard_terminus = FALSE,
            cluster_labels=FALSE)
  {
    if (!is.null(clusters))
      x$metaclusters <- clusters
    if (is.null(x$metaclusters)) {
      message("Computing louvain metaclusters...\n")
      if (is.null(x$dsym)) {
        d <- KofRawN(x$KNN, K)
        d <- knn.raw2adj(d)
        dsym <- knn.spadj2sym(knn.adj2spadj(d))
      }
      else
        dsym <- x$dsym
      gU <- graph_from_adjacency_matrix(dsym, weighted = TRUE,
                                        mode = "undirected")
      x$metaclusters <- cluster_louvain(gU)$membership
      message("~Done!\n")
      rm(gU)
    }
    if (is.null(x$dsym)) {
      d <- KofRawN(x$KNN, K)
      d <- knn.raw2adj(d)
      dsym <- knn.spadj2sym(knn.adj2spadj(d))
    }
    else
      dsym <- x$dsym
    e.list <-
      cbind(x$walks[[origin_name]]$v[-c(x$walks[[origin_name]]$starts[-1] -
                                          1, length(x$walks[[origin_name]]$v))], x$walks[[origin_name]]$v[-x$walks[[origin_name]]$starts])
    e.list <- rbind(e.list, dim(dsym))
    gD <- graph_from_edgelist(e.list, directed = TRUE)
    gD <- contract(gD, x$metaclusters)
    agDc <- summary(as_adjacency_matrix(gD))
    sel <- which(agDc$i == agDc$j)
    agDc <- agDc[-sel,]
    clus <- sort(unique(x$metaclusters))
    agDc <- sparseMatrix(
      i = agDc$i,
      j = agDc$j,
      x = agDc$x,
      dims = rep(length(clus), 2)
    )
    g_layout <- NULL
    ins <- Matrix::rowSums(Matrix::t(agDc))
    outs <- Matrix::rowSums(agDc)
    ratio <- outs / ins
    ratio[ins == 0] <- 1
    if (hard_terminus)
      agDc[ratio < 1, ] <- 0
    S <- x$metaclusters[x$origin[[origin_name]]]
    for (i in clus)
      g_layout <-
      rbind(g_layout, colMeans(x$layout[[layout]][which(x$metaclusters ==
                                                          i),]))
    for (i in 1:dim(agDc)[1]) {
      for (j in 1:dim(agDc)[1]) {
        if (agDc[i, j] < agDc[j, i])
          agDc[i, j] <- 0
        else
          agDc[j, i] <- 0
      }
    }
    if (scale == "max") {
      scale_factor <- (ins + outs) / 2
      scale_factor <- rep(max(scale_factor), length(scale_factor))
    }
    else if (scale == "adaptive") {
      scale_factor <- (ins + outs) / 2
      scale_factor[which(ratio < 1)] <- ins[which(ratio < 1)]
      S <- x$metaclusters[x$origin[[origin_name]]]
      scale_factor[S] <- outs[S]
    }
    else {
      scale_factor <- (ins + outs) / 2
    }
    aa1 <- Diagonal(x = scale_factor ^ -1) %*% agDc
    g11 <- graph_from_adjacency_matrix(aa1, weighted = TRUE,
                                       mode = "directed")
    E(g11)$width <- E(g11)$weight ^ 1.2 * 19
    E(g11)$curved = TRUE
    E(g11)$arrow.size <- 2.7 * arrow.sizefactor
    E(g11)$arrow.width <- 1 * arrow.sizefactor + 1
    E(g11)$arrow.mode <- rep(2, length(E(g11)))
    E(g11)$arrow.mode[E(g11)$width < quantile(E(g11)$width, qq)] <-
      0
    E(g11)$width[E(g11)$width > quantile(E(g11)$width, qq)] <-
      E(g11)$width[E(g11)$width > quantile(E(g11)$width, qq)] + 7
    V(g11)$label.cex = 4
    V(g11)$label <- 1:length(clus)
    if (!cluster_labels) V(g11)$label <- rep("", length(clus))
    V(g11)$label[ratio < 1] <-
      paste("T", V(g11)$label[which(ratio < 1)])
    V(g11)$label[S] <- paste("O", V(g11)$label[S])
    pieD <- list()[1:length(clus)]
    for (i in clus) {
      if (!is.null(notplot)) {
        sel <- which(!(x$labels[[labels]][x$metaclusters ==
                                            i] %in% notplot))
        pieD[[i]] <-
          as.vector(table(x$labels[[labels]][x$metaclusters ==
                                               i][sel]))
        if (length(pieD[[i]]) == 0)
          pieD[[i]] <- length(unique(x$labels[[labels]])) +
          1
      }
      else {
        pieD[[i]] <- as.vector(table(x$labels[[labels]][x$metaclusters ==
                                                          i]))
      }
    }
    cp <-
      c(
        RColorBrewer::brewer.pal(8, "Dark2"),
        RColorBrewer::brewer.pal(12,
                                 "Paired")[-11],
        RColorBrewer::brewer.pal(9, "Set1")[-6],
        RColorBrewer::brewer.pal(8, "Accent")[5:8]
      )
    colors <- list(cp)
    frame.color <- rep("black", length(clus))
    frame.color[S] <- "gold"
    frame.color[ratio < 1] <- "darkred"
    vertex.size <- rep(7, length(clus))
    vertex.size[ratio < 1] <- 12
    vertex.size[S] <- 12
    colors <- rep(colors, length(clus))
    lsize = 0.4
    if (!is.null(png)) {
      png(png, 2000, 2000)
      lsize = 2.3
    }
    igraph::plot.igraph(
      g11,
      layout = g_layout,
      main = "infered connectome",
      vertex.size = vertex.size,
      vertex.shape = "pie",
      vertex.pie = pieD,
      vertex.pie.color = colors,
      vertex.frame.color = frame.color,
      vertex.pie.border = frame.color,
      vertex.frame.width = 10
    )
    
    ss <- which(levels(x$labels[[labels]]) != notplot)
    if (length(ss) == 0)
      ss <- 1:length(x$labels[[labels]])
    legend(
      legend.position,
      legend = levels(x$labels[[labels]])[ss],
      col = colors[[1]][ss],
      pch = 19,
      cex = lcex
    )
    if (!is.null(png))
      dev.off()
    return(invisible(x))
  }