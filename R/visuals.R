# plot.walks <- function(
#   X = proj.clusters$layout,
#   w = walks.clusters,
#   sel = walks.selection,
#   col = 1:length(walks$starts),
#   ...
# ) {
#   j <- 0
#   for (i in sel) {
#     j <- j + 1
#     pts <- X[select_paths_points(w, i), ]
#     lines(pts, lty = 1, col = col[j], ...)
#     points(pts[, 1], pts[, 2], pch = 20, cex = 1.5, col = col[j])
#   }
# }
# 
# plot.walks.list <- function(
#   X = proj.C$layout,
#   pathways = pathways,
#   col = 1:length(unique(unlist(pathways$classes))),
#   ...
# ) {
#   j <- 0
#   for (i in 1:length(pathways$walks)) {
#     j <- j + 1
#     pts <- X[pathways$walks[[i]], ]
#     lines(pts, lty = 1, col = col[pathways$classes[[i]]], ...)
#     points(pts[, 1], pts[, 2], pch = 20, cex = 1.5, col = col[pathways$classes[[i]]])
#   }
# }
# 
# plot.consensus_paths <- function(
#   cons = consensus_paths,
#   col = 1:length(cons),
#   ...
# ) {
#   j <- 0
#   for (i in 1:length(cons)) {
#     j <- j + 1
#     pts <- cons[[i]]
#     lines(pts, lty = 1, col = col[i], ...)
#     points(pts[, 1], pts[, 2], pch = 20, cex = 1.5, col = col[i])
#   }
# }

#' Visualise pathways progression: single parameter by segments
#'
#' \code{pathways.visual.parameters.segments} generates a \code{ggplot2} visualisation of differential expression
#'  for a single parameter with pseudotime progression, grouped by pathways coming from a single class. Pseudotime
#'  progression is divided into segments for which mean expression per pathway is computed and plotted.
#' @param paths a list of pathways. This corresponds to the \code{walks} slot of a pathways list object, as generated
#'  by \code{create_pathways_object}.
#' @param classes a list of classes corresponding to pathways. This corresponds to the \code{classes} slot of a pathways list object, as generated
#'  by \code{create_pathways_object}.
#' @param which.class number of the pathway class used for generating the plot.
#' @param coordinates a coordinate matrix to which the input pathways refer. This is a matrix of expression values, column correspond to
#'  parameters (markers) and rows correspond to single events.
#' @param pseudotime a pseudotime value vector, corresponding to vertices of \code{coordinates}.
#' @param which.param name of the parameter (from \code{coordinates} column names) used for generating the plot.
#' @param breaks if not null, this overrides \code{n.part} and \code{exp.part}: create pseudotime partitions by dividing the
#'  \code{pseudotime} vector by its values (instead of including a fixed number of cells in each segment). \codes{breaks} can
#'  consist either of a single values (number of uniform partitions) or of a vector of breaks (cut-off values) between 0 and 1.
#' @param n.part number of partitions per pseudotime. Used which \code{breaks} is null. A fixed number of cells, picked from a
#'  series sorted by pseudotime values, is included in each partition.
#' @param exp.part used with \code{n.part} if \code{breaks} is null. Non-negative exponent used in computing the series of cut-off values for
#'  pseudotime vector partitioning. If equal to 1, partitioning is uniform (each parition contains an equal number of cells). If
#'  greater than one, partition size shrinks with increasing pseudotime (recommended, useful if many mature cells/terminal events are present in sample).
#'  If between 0 and 1, partition size expands with increasing pseudotime.
#' @param expression.transform transformation to apply on expression values (entries of \code{coordinates}). Accepts \code{"asinh"}.
#'  Defaults to null (no transformation).
#' @param asinh.denominator if \code{expression.transform} is \code{"asinh"}, which factor \code{F} to use in \code{coordinates <- asinh(coordinates/F)}.
#'  \code{F = 5} recommended for CyTOF data, approximately \code{F = 120} for flow data.
#'
#' @return \code{pathways.visual.parameters.segments} returns a \code{ggplot} object.
#'
#' @export
pathways.visual.parameters.segments <- function(
  paths,
  classes,
  which.class,
  coordinates,
  pseudotime,
  which.param,
  breaks = NULL,
  n.part = 10,
  exp.part = 1,
  expression.transform = NULL,
  asinh.denominator = 5
) {
  require(ggplot2)
  
  pseudotime <- pseudotime/max(pseudotime)
  
  classes.inds <- classes %in% which.class
  params.inds <- colnames(coordinates) %in% which.param
  
  coordinates <- coordinates[, params.inds]
  
  if (!is.null(expression.transform)) {
    if (expression.transform == "asinh") {
      coordinates <- asinh(coordinates/asinh.denominator)
    }
  }
  
  paths <- paths[classes.inds]
  progress <- lapply(paths, function(pts) {
    pseudotime[pts] # pseudotime progress per path
  })
  
  pp <- unique(unlist(progress))
  
  b <- NULL
  N <- NULL
  
  if (is.null(breaks)) {
    p <- (seq(0, 1, 1/(n.part)))^exp.part
    b <- as.numeric(quantile(pp, p))
    N <- n.part
  } else {
    if (length(breaks) > 1) {
      b <- breaks
      N <- max(breaks)
    } else {
      b <- sapply(0:breaks, function(i) {
        1/breaks * i
      })
      N <- length(b) - 1
    }
  }
  
  categs <- as.numeric(cut(unlist(progress), breaks = b, include.lowest = TRUE))
  
  stats <- lapply(1:N, function(i) {
    inds <- categs == i # pick points on paths by pseudotime increment
    
    if (!any(inds)) {
      return(NULL)
    }
    
    pts <- unlist(paths)[inds]
    
    which.paths <- unlist(lapply(1:length(paths), function(j) {
      rep(j, length(paths[[j]]))
    }))[inds]
    
    pts <- unlist(paths)[unlist(inds)]
    
    means <- lapply(sort(unique(which.paths)), function(j) {
      mean(coordinates[pts[which.paths == j]])
    })
    
    v <- cbind(rep(i, length(means)), sort(unique(which.paths)), unlist(means))
    if (is.null(dim(v))) {
      names(v) <- c("segment", "path", "expression")
    } else {
      colnames(v) <- c("segment", "path", "expression")
    }
    v <- as.data.frame(v)
    v
  })
  stats <- do.call(rbind, stats)
  
  stats$path <- as.factor(stats$path)
  stats$segment <- as.numeric(stats$segment)
  stats$expression <- as.numeric(stats$expression)
  
  ggplot(data = stats, aes(x = segment, y = expression, group = path, color = path)) + 
    geom_line(linetype = "dashed") +
    geom_point() +
    ggtitle(paste(which.param, " expression per path of class ", which.class, ": means and standard deviations", sep = "")) +
    labs(subtitle = paste(
      "Segments were generated based on pseudotime values. ",
      ifelse(!is.null(expression.transform),
             paste("\n", expression.transform, " transformation was applied to expression values",
                   ifelse(expression.transform == "asinh",
                          paste(", with denominator ", asinh.denominator, ".", sep = ""),
                          "."
                   ), sep = ""),
             ""
      ), sep = "")) +
    theme_minimal()
}

#' Visualise pathways progression: single parameter using all pathway points
#'
#' \code{pathways.visual.parameters.all_data} generates a \code{ggplot2} visualisation of differential expression
#'  for a single parameter with pseudotime progression, showing all pathways coming from a single class.
#' @param paths a list of pathways. This corresponds to the \code{walks} slot of a pathways list object, as generated
#'  by \code{create_pathways_object}.
#' @param classes a list of classes corresponding to pathways. This corresponds to the \code{classes} slot of a pathways list object, as generated
#'  by \code{create_pathways_object}.
#' @param which.class number of the pathway class used for generating the plot.
#' @param coordinates a coordinate matrix to which the input pathways refer. This is a matrix of expression values, column correspond to
#'  parameters (markers) and rows correspond to single events.
#' @param pseudotime a pseudotime value vector, corresponding to vertices of \code{coordinates}.
#' @param which.param name of the parameter (from \code{coordinates} column names) used for generating the plot.
#' @param expression.transform transformation to apply on expression values (entries of \code{coordinates}). Accepts \code{"asinh"}.
#'  Defaults to null (no transformation).
#' @param asinh.denominator if \code{expression.transform} is \code{"asinh"}, which factor \code{F} to use in \code{coordinates <- asinh(coordinates/F)}.
#'  \code{F = 5} recommended for CyTOF data, approximately \code{F = 120} for flow data.
#'
#' @return \code{pathways.visual.parameters.segments} returns a \code{ggplot} object.
#'
#' @export
pathways.visual.parameters.all_data <- function(
  paths,
  classes,
  which.class,
  coordinates,
  pseudotime,
  which.param,
  expression.transform = NULL,
  asinh.denominator = 5
) {
  require(ggplot2)
  
  pseudotime <- pseudotime/max(pseudotime)
  
  classes.inds <- classes %in% which.class
  params.inds <- colnames(coordinates) %in% which.param
  
  coordinates <- coordinates[, params.inds]
  
  if (!is.null(expression.transform)) {
    if (expression.transform == "asinh") {
      coordinates <- asinh(coordinates/asinh.denominator)
    }
  }
  
  paths <- paths[classes.inds]
  progress <- lapply(paths, function(pts) {
    pseudotime[pts] # pseudotime progress per path
  })
  which.path <- lapply(1:length(progress), function(i) {
    rep(i, length(progress[[i]]))
  })
  
  m <- max(unlist(progress))
  progress <- lapply(progress, function(p) {
    p <- p/m
  })
  
  expression <- lapply(paths, function(p) {
    coordinates[p]
  })
  
  d <- as.data.frame(cbind(unlist(which.path), unlist(progress), unlist(expression))); colnames(d) <- c("path", "pseudotime", "expression")
  d$path <- as.factor(d$path)
  
  g <- ggplot(data = d, aes(x = pseudotime, y = expression, group = path, color = path)) + 
    geom_line() +
    ggtitle(paste(which.param, "expression per path of class ", which.class, ": individual events", sep = "")) +
    labs(subtitle = paste("Differential expression across pseudotime is shown, grouped by paths. ",
                          ifelse(!is.null(expression.transform),
                                 paste("\n", expression.transform, " transformation was applied to expression values",
                                       ifelse(expression.transform == "asinh",
                                              paste(", with denominator ", asinh.denominator, ".", sep = ""),
                                              "."
                                       ), sep = ""),
                                 ""
                          )
    )) +
    theme_minimal()
  
  g
}
#' Visualise pathways progression: cell population labels by segment
#'
#' \code{pathways.visual.parameters.segments} generates a \code{ggplot2} visualisation of differential representation
#'  of labelled cell populations with pseudotime progression, summing out pathways from a single class over pseudotime segments.
#' @param paths a list of pathways. This corresponds to the \code{walks} slot of a pathways list object, as generated
#'  by \code{create_pathways_object}.
#' @param classes a list of classes corresponding to pathways. This corresponds to the \code{classes} slot of a pathways list object, as generated
#'  by \code{create_pathways_object}.
#' @param which.class number of the pathway class used for generating the plot.
#' @param labels cell population labels (eg. gates) corresponding to single events (to which \code{paths} refer).
#' @param pseudotime a pseudotime value vector, corresponding to vertices of \code{coordinates}.
#' @param breaks if not null, this overrides \code{n.part} and \code{exp.part}: create pseudotime partitions by dividing the
#'  \code{pseudotime} vector by its values (instead of including a fixed number of cells in each segment). \codes{breaks} can
#'  consist either of a single values (number of uniform partitions) or of a vector of breaks (cut-off values) between 0 and 1.
#' @param n.part number of partitions per pseudotime. Used which \code{breaks} is null. A fixed number of cells, picked from a
#'  series sorted by pseudotime values, is included in each partition.
#' @param exp.part used with \code{n.part} if \code{breaks} is null. Non-negative exponent used in computing the series of cut-off values for
#'  pseudotime vector partitioning. If equal to 1, partitioning is uniform (each parition contains an equal number of cells). If
#'  greater than one, partition size shrinks with increasing pseudotime (recommended, useful if many mature cells/terminal events are present in sample).
#'  If between 0 and 1, partition size expands with increasing pseudotime.
#' @param proportional Boolean. If \code{TRUE}, show relative representations of cell populations. If \code{FALSE}, show absolute cell counts.
#'  Default to \code{TRUE}.
#'
#' @return \code{pathways.visual.labels.segments} returns a \code{ggplot} object.
#'
#' @export
pathways.visual.labels.segments <- function(
  paths,
  classes,
  which.class,
  labels,
  pseudotime,
  breaks = NULL,
  n.part = 20,
  exp.part = 1,
  proportional = TRUE
) {
  require(ggplot2)
  
  pseudotime <- pseudotime/max(pseudotime)
  
  classes.inds <- classes %in% which.class
  
  paths <- paths[classes.inds]
  progress <- lapply(paths, function(pts) {
    pseudotime[pts] # pseudotime progress per path
  })
  
  gates <- lapply(paths, function(pts) {
    labels[pts]
  })
  gates.u <- unique(labels)
  
  pp <- unique(unlist(progress))
  
  b <- NULL
  N <- NULL
  
  if (is.null(breaks)) {
    p <- (seq(0, 1, 1/(n.part)))^exp.part
    b <- as.numeric(quantile(pp, p))
    N <- n.part
  } else {
    if (length(breaks) > 1) {
      b <- breaks
      N <- max(breaks)
    } else {
      b <- sapply(0:breaks, function(i) {
        1/breaks * i
      })
      N <- length(b) - 1
    }
  }
  
  categs <- as.numeric(cut(unlist(progress), breaks = b, include.lowest = TRUE))
  
  stats <- lapply(1:N, function(i) {
    
    inds <- categs == i # pick points on paths by pseudotime increment
    
    if (!any(inds)) {
      return(NULL)
    }
    
    g <- (as.data.frame(table(unlist(gates)[inds]))); colnames(g) <- c("label", "count")
    g <- cbind(i, g); colnames(g)[1] <- "segment"
    
    missing <- gates.u[!gates.u %in% g$label]
    if (length(missing) > 0) {
      missing <- as.data.frame(do.call(rbind, lapply(missing, function(m) c(i, m, 0)))); colnames(missing) <- c("segment", "label", "count")
      g <- rbind(g, missing)
    }
    
    g$segment <- as.numeric(g$segment)
    g$label <- as.factor(g$label)
    g$count <- as.numeric(g$count)
    
    if (proportional) {
      g$count <- g$count / sum(g$count)
    }
    
    g
  })
  stats <- do.call(rbind, stats)
  
  stats$label <- factor(stats$label, levels = sort(levels(l)))
  
  g <- ggplot(stats, aes(x = segment, y = count, fill = label)) + 
    geom_area(alpha = .9, colour = "black", size = .2) +
    labs(fill = "Labels", alpha = .9) +
    ggtitle(paste("Labels in path of class ", which.class, ": pseudotime segments", sep = "")) +
    labs(subtitle = paste("Differential representation of populations across pseudotime is shown."
    )) +
    theme_minimal()
  
  if (proportional) {
    g <- g + ylab("representation")
  }
  
  g
}
