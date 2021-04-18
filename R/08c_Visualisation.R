.draw_placeholder <- function(
  picture = 'tree'
) {
  if (picture == 'tree') {
    j <- jpeg::readJPEG(system.file('tree.jpg', package = 'tviblindi'), native = TRUE)
  }
  if (picture == 'petal') {
    j <- jpeg::readJPEG(system.file('petal.jpg', package = 'tviblindi'), native = TRUE)
  }
  plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
  rasterImage(j, 0, 0, 1, 1)
}

.plot_trajectories <- function(
  X,
  walks,
  walk_idcs.A,
  walk_idcs.B,
  flip_colours,
  pseudotime_highlight_bounds,
  pseudotime,
  highlight_in_background,
  selected_trajectory_points = NULL,
  pointsize                  = .05,
  ...
) {
  ## Plot trajectories over a 2D layout
  col1 <- c(34, 87, 201, 255)
  col2 <- c(194, 45, 55, 255)
  plot(scattermore::scattermore(X, rgba = c(200, 200, 200, 150), cex = pointsize))
  
  j      <- 0
  sel1   <- walks[walk_idcs.A]
  sel2   <- walks[walk_idcs.B]
  
  if (!is.null(pseudotime_highlight_bounds) && highlight_in_background) {
    p <- pseudotime$res / max(pseudotime$res)
    a <- which(p >= pseudotime_highlight_bounds[1])
    b <- which(p <= pseudotime_highlight_bounds[2])
    idcs.highlight <- intersect(a, b)
    idcs.highlight <- idcs.highlight[idcs.highlight %in% unlist(c(sel1, sel2))]
    
    if (length(idcs.highlight) > 0) {
      
      pts <- X[idcs.highlight, , drop = FALSE]
      
      if (!is.null(selected_trajectory_points)) pts<-X[selected_trajectory_points, ,drop=FALSE]
      
      plot(scattermore::scattermore(pts, rgba = c(0, 153, 31, 255), xlim = c(0, 1), ylim = c(0, 1), cex = pointsize + 1), add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
    }
  }
  
  if (flip_colours) {
    tmp <- sel1
    sel1 <- sel2
    sel2 <- tmp
    tmp <- col1
    col1 <- col2
    col2 <- tmp
  }
  
  if (length(sel1) > 0) {
    pts1 <- vector(mode = 'list', length = length(sel1))
    for (idx in 1:length(sel1)) {
      s   <- sel1[[idx]]
      j   <- j + 1
      pts <- X[s, ]
      
      pts1[[idx]] <- do.call(rbind, interpolate_trajectories(lapply(1:nrow(pts), function(x) pts[x, ])))
    }
    pts1 <- do.call(rbind, pts1)
    plot(scattermore::scattermore(pts1, rgba = col1, xlim = c(0, 1), ylim = c(0, 1)), add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
  }
  
  if (length(sel2) > 0) {
    pts2 <- vector(mode = 'list', length = length(sel2))
    for (idx in 1:length(sel2)) {
      s   <- sel2[[idx]]
      j   <- j + 1
      pts <- X[s, ]
      
      pts2[[idx]] <- do.call(rbind, interpolate_trajectories(lapply(1:nrow(pts), function(x) pts[x, ])))
    }
    pts2 <- do.call(rbind, pts2)
    plot(scattermore::scattermore(pts2, rgba = col2, xlim = c(0, 1), ylim = c(0, 1)), add = TRUE)
  }
  
  if (!is.null(pseudotime_highlight_bounds) && !highlight_in_background) {
    p <- pseudotime$res / max(pseudotime$res)
    a <- which(p >= pseudotime_highlight_bounds[1])
    b <- which(p <= pseudotime_highlight_bounds[2])
    idcs.highlight <- intersect(a, b)
    idcs.highlight <- idcs.highlight[idcs.highlight %in% unlist(c(sel1, sel2))]
    
    if (length(idcs.highlight) > 0) {
      pts <- X[idcs.highlight, , drop = FALSE]
      if (!is.null(selected_trajectory_points)) pts<-X[selected_trajectory_points, ,drop=FALSE]
      plot(scattermore::scattermore(pts, rgba = c(0, 153, 31, 255), xlim = c(0, 1), ylim = c(0, 1), cex = pointsize + 1), add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
    }
  }
}

.plot_trajectories_full <- function(
  X,
  walks,
  walk_idcs.A,
  walk_idcs.B,
  flip_colours,
  pseudotime_highlight_bounds,
  pseudotime,
  highlight_in_background,
  pointsize = .08,
  ...
) {
  cols <- c('blue', 'darkred')
  j      <- 0
  sel1   <- walks[walk_idcs.A]
  sel2   <- walks[walk_idcs.B]
  
  plot(X, pch = 20, cex = pointsize, col = 'darkgrey', xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = '', ylab = '')
  
  if (!is.null(pseudotime_highlight_bounds) && highlight_in_background) {
    p <- pseudotime$res / max(pseudotime$res)
    a <- which(p >= pseudotime_highlight_bounds[1])
    b <- which(p <= pseudotime_highlight_bounds[2])
    idcs.highlight <- intersect(a, b)
    idcs.highlight <- idcs.highlight[idcs.highlight %in% unlist(c(sel1, sel2))]
    
    if (length(idcs.highlight) > 0) {
      pts <- X[idcs.highlight, , drop = FALSE]
      points(pts, col = 'lightgreen', pch = 20, cex = pointsize + 1)
    }
  }
  
  if (flip_colours) {
    tmp <- sel1
    sel1 <- sel2
    sel2 <- tmp
    cols <- cols[c(2, 1)]
  }
  i <- c(walk_idcs.A, walk_idcs.B)
  a.idcs <- if (is.null(i)) { (1:length(walks)) } else { (1:length(walks))[-i] }
  
  rest   <- if (length(a.idcs) > 0) { walks[a.idcs] } else { NULL }
  for (i in rest) {
    j   <- j + 1
    pts <- X[i, ]
    lines(pts, lty = 1, col = scales::alpha('darkgrey', 0.02), ...)
    points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = scales::alpha('darkgrey', 0.02))
  }
  if (length(sel1) > 0) {
    for (i in sel1) {
      j   <- j + 1
      pts <- X[i, ]
      lines(pts, lty = 1, col = scales::alpha(cols[1], 0.3), ...)
      points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = scales::alpha(cols[1], 0.3))
    }
  }
  if (length(sel2) > 0) {
    for (i in sel2) {
      j   <- j + 1
      pts <- X[i, ]
      lines(pts, lty = 1, col = scales::alpha(cols[2], 0.3), ...)
      points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = scales::alpha(cols[2], 0.3))
    }
  }
  if (!is.null(pseudotime_highlight_bounds) && !highlight_in_background) {
    p <- pseudotime$res / max(pseudotime$res)
    a <- which(p >= pseudotime_highlight_bounds[1])
    b <- which(p <= pseudotime_highlight_bounds[2])
    idcs.highlight <- intersect(a, b)
    idcs.highlight <- idcs.highlight[idcs.highlight %in% unlist(c(sel1, sel2))]
    
    if (length(idcs.highlight) > 0) {
      pts <- X[idcs.highlight, , drop = FALSE]
      points(pts, col = 'lightgreen', pch = 20, cex = pointsize + 1)
    }
  }
}

## Functions: plot expression markers across pseudotime for selected walks
.plot_tracked_markers <- function(
  walks,
  walk_idcs,
  tv,
  pseudotime,
  markers,
  breaks          = NULL,
  n.part          = 10,
  exp.part        = 1,
  large_base_size = FALSE,
  grey            = FALSE
) {
  if (length(markers) == 1) {
    return(do.call(.plot_tracked_markers_single, as.list(environment())))
  } else {
    return(do.call(.plot_tracked_markers_multiple, as.list(environment())))
  }
}

.plot_tracked_markers_single <- function(
  walks,
  walk_idcs,
  tv,
  pseudotime,
  markers,
  breaks = NULL,
  n.part = 10,
  exp.part = 1,
  show_legend = FALSE,
  large_base_size = FALSE,
  grey = FALSE
) {
  ## Prevent 'no visible binding for global variable' build warnings
  segment <- marker <- walk <- count <- NULL 
  
  ## Scale pseudotime
  pseudotime <- pseudotime$res
  pseudotime <- pseudotime/max(pseudotime)
  
  ## Pick relevant walks
  walks      <- walks[walk_idcs]
  
  ## Get pseudotimes per walk point for each walk
  progress   <- lapply(walks, function(pts) pseudotime[pts])
  pp         <- unique(unlist(progress))
  
  b <- NULL
  N <- NULL
  
  if (is.null(breaks)) {
    p <- (seq(0, 1, 1 / (n.part))) ^ exp.part # re-scale timeline of pseudotime segments
    b <- as.numeric(stats::quantile(pp, p))
    N <- n.part
  } else {
    if (length(breaks) > 1) {
      b <- breaks
      N <- max(breaks)
    } else {
      b <- sapply(0:breaks, function(i) 1 / breaks * i)
      N <- length(b) - 1
    }
  }
  
  categs <- as.numeric(cut(unlist(progress), breaks = b, include.lowest = TRUE))
  # for each point on walk: which segment does it fall into?
  
  pseudotime_bounds <- sort(c(unlist(progress)[which(!duplicated(categs))], 1))
  coords <- tv$data[, markers]
  
  stats  <- lapply(1:N, function(i) {
    inds <- categs == i # pick points on paths by pseudotime increment
    if (!any(inds)) return(NULL)
    
    which.walks <- unlist(lapply(1:length(walks), function(j) rep(j, length(walks[[j]]))))[inds]
    pts <- unlist(walks)[unlist(inds)]
    means <- lapply(sort(unique(which.walks)), function(j) mean(coords[pts[which.walks == j]]))
    
    inds_char <- lapply(sort(unique(which.walks)), function(j) paste(pts[which.walks == j], collapse =","))
    
    v <- cbind(rep(i, length(means)), sort(unique(which.walks)), unlist(means))
    if (is.null(dim(v))) { names(v) <- c('segment', 'walk', 'expression') } else { colnames(v) <- c('segment', 'walk', 'expression') }
    as.data.frame(v)
  })
  
  inds_char  <- lapply(1:N, function(i) {
    inds <- categs == i # pick points on paths by pseudotime increment
    if (!any(inds)) return(NULL)
    which.walks <- unlist(lapply(1:length(walks), function(j) rep(j, length(walks[[j]]))))[inds]
    pts         <- unlist(walks)[unlist(inds)]
    inds_char  <- lapply(sort(unique(which.walks)), function(j) paste(pts[which.walks == j],collapse =","))
    return(matrix(unlist(inds_char),ncol=1))
  })
  
  inds_char<-do.call(rbind,inds_char)
  stats            <- do.call(rbind, stats)
  stats$walk       <- as.factor(stats$walk)
  stats$segment    <- as.numeric(stats$segment)
  stats$expression <- as.numeric(stats$expression)
  
  g <- ggplot(data = stats, aes(x = segment, y = expression, group = walk, color = walk)) +
    geom_line(linetype = 'dashed') + geom_point() +
    ggtitle(paste0(markers, ' expression per walk: means per segment')) +
    labs(subtitle = 'Segmented by pseudotime values.') +
    theme(text          = element_text(size = if (large_base_size) { 24 } else { 11 }),
          axis.text.x   = element_blank(),
          axis.ticks    = element_blank(),
          plot.title    = element_text(size = if (large_base_size) { 20 } else { 18 }),
          plot.subtitle = element_text(size = if (large_base_size) { 19 } else { 16 }),
          legend.title  = element_text(size = if (large_base_size) { 24 } else { 16 }),
          legend.text   = element_text(size = if (large_base_size) { 22 } else { 12 }))
  if (grey) {
    g <- g + theme(panel.background = element_rect(fill = '#f2f2f2',
                                                   colour = '#f2f2f2',
                                                   size = 0.5, linetype = 'solid'))
  } else {
    g <- g + theme_minimal()
  }
  if (!show_legend) {
    g <- g + theme(legend.position = 'none')
  }
  list(plot              = g,
       stats             = cbind(stats,inds_char=inds_char),
       pseudotime_bounds = pseudotime_bounds)
}

.plot_tracked_markers_multiple <- function(
  walks,
  walk_idcs,
  tv,
  pseudotime,
  markers,
  breaks          = NULL,
  n.part          = 10,
  exp.part        = 1,
  show_legend     = TRUE,
  large_base_size = FALSE,
  grey            = FALSE
) {
  ## Prevent 'no visible binding for global variable' build warnings
  segment <- marker <- walk <- count <- NULL 
  
  pseudotime <- pseudotime$res
  pseudotime <- pseudotime/max(pseudotime)
  
  walks    <- walks[walk_idcs]
  progress <- lapply(walks, function(pts) pseudotime[pts])
  pp       <- unique(unlist(progress))
  
  b <- NULL
  N <- NULL
  
  if (is.null(breaks)) {
    p <- (seq(0, 1, 1 / (n.part))) ^ exp.part
    b <- as.numeric(stats::quantile(pp, p))
    N <- n.part
  } else {
    if (length(breaks) > 1) {
      b <- breaks
      N <- max(breaks)
    } else {
      b <- sapply(0:breaks, function(i) 1 / breaks * i)
      N <- length(b) - 1
    }
  }
  
  categs <- as.numeric(cut(unlist(progress), breaks = b, include.lowest = TRUE))
  
  pseudotime_bounds <- sort(c(unlist(progress)[which(!duplicated(categs))], 1))
  
  coords <- tv$data[, markers]
  
  stats  <- lapply(markers, function(m) {
    lapply(1:N, function(i) {
      
      inds <- categs == i # pick points on paths by pseudotime increment
      if (!any(inds)) return(NULL)
      
      which.walks <- unlist(lapply(1:length(walks), function(j) rep(j, length(walks[[j]]))))[inds]
      pts         <- unlist(walks)[unlist(inds)]
      means       <- sapply(sort(unique(which.walks)), function(j) mean(coords[pts[which.walks == j], m])) %>% mean
      
      
      v <- cbind(rep(m, length(means)), rep(i, length(means)), unlist(means))
      if (is.null(dim(v))) { names(v) <- c('marker', 'segment', 'expression') } else { colnames(v) <- c('marker', 'segment', 'expression') }
      as.data.frame(v)
    })
  })
  stats            <- unlist(stats, recursive = FALSE)
  stats            <- do.call(rbind, stats)
  stats$marker     <- as.factor(stats$marker)
  stats$segment    <- as.numeric(stats$segment)
  stats$expression <- as.numeric(as.character(stats$expression))
  
  g <- ggplot(data = stats, aes(x = segment, y = expression, group = marker, color = marker)) +
    geom_line() +
    geom_point() +
    ggtitle(paste0('Multiple markers expression')) +
    labs(subtitle = 'Segmented by pseudotime values. ') +
    theme(
      text          = element_text(size = if (large_base_size) { 24 } else { 11 }),
      axis.text.x   = element_blank(),
      axis.ticks    = element_blank(),
      plot.title    = element_text(size = if (large_base_size) { 20 } else { 18 }),
      plot.subtitle = element_text(size = if (large_base_size) { 19 } else { 16 }),
      legend.title  = element_text(size = if (large_base_size) { 24 } else { 16 }),
      legend.text   = element_text(size = if (large_base_size) { 22 } else { 12 })
    )
  if (grey) {
    g <- g + theme(panel.background = element_rect(
      fill     = '#f2f2f2',
      colour   = '#f2f2f2',
      size     = 0.5,
      linetype = 'solid')
    )
  } else {
    g <- g + theme_minimal()
  }
  list(
    plot              = g,
    stats             = cbind(stats,inds_char="NULL"),
    pseudotime_bounds = pseudotime_bounds
  )
}

.plot_tracked_populations <- function(
  walks,
  walk_idcs,
  tv,
  labels_name,
  pseudotime,
  populations,
  breaks          = NULL,
  n.part          = 20,
  exp.part        = 1,
  large_base_size = FALSE,
  log2_transform  = NULL,
  grey            = FALSE
) {
  ## Prevent 'no visible binding for global variable' build warnings
  segment <- walk <- count <- population <- NULL 
  
  pseudotime <- pseudotime$res
  pseudotime <- pseudotime / max(pseudotime)
  
  walks    <- walks[walk_idcs]
  progress <- lapply(walks, function(pts) pseudotime[pts])
  pp       <- unique(unlist(progress))
  
  b <- NULL
  N <- NULL
  
  if (is.null(breaks)) {
    p <- (seq(0, 1, 1 / (n.part))) ^ exp.part
    b <- as.numeric(stats::quantile(pp, p))
    N <- n.part
  } else {
    if (length(breaks) > 1) {
      b <- breaks
      N <- max(breaks)
    } else {
      b <- sapply(0:breaks, function(i) 1 / breaks * i)
      N <- length(b) - 1
    }
  }
  
  categs <- as.numeric(cut(unlist(progress), breaks = b, include.lowest = TRUE))
  
  pseudotime_bounds <- sort(unlist(progress)[which(!duplicated(categs))])
  
  stats  <- lapply(populations, function(p) {
    lapply(1:N, function(i) {
      
      inds <- categs == i # pick points on paths by pseudotime increment
      if (!any(inds)) return(NULL)
      
      which.walks <- unlist(lapply(1:length(walks), function(j) rep(j, length(walks[[j]]))))[inds]
      pts         <- unique(unlist(walks)[unlist(inds)])
      counts      <- sum(tv$labels[[labels_name]][pts] == p)
      if (log2_transform) {
        counts <- log2(counts)
      }
      
      v <- cbind(rep(p, length(counts)), rep(i, length(counts)), unlist(counts))
      if (is.null(dim(v))) { names(v) <- c('population', 'segment', 'count') } else { colnames(v) <- c('population', 'segment', 'count') }
      as.data.frame(v)
    })
  })
  stats            <- unlist(stats, recursive = FALSE)
  stats            <- do.call(rbind, stats)
  stats$population <- as.factor(stats$population)
  stats$segment    <- as.numeric(stats$segment)
  stats$count      <- as.numeric(as.character(stats$count))
  
  g <- ggplot(data = stats, aes(x = segment, y = count, group = population, color = population)) +
    geom_line() +
    geom_point() +
    ggtitle(paste0('Annotated populations composition progression')) +
    labs(subtitle = 'Segmented by pseudotime values. ') +
    theme(text          = element_text(size = if (large_base_size) { 24 } else { 11 }),
          axis.text.x   = element_blank(),
          axis.ticks    = element_blank(),
          plot.title    = element_text(size = if (large_base_size) { 20 } else { 18 }),
          plot.subtitle = element_text(size = if (large_base_size) { 19 } else { 16 }),
          legend.title  = element_text(size = if (large_base_size) { 24 } else { 16 }),
          legend.text   = element_text(size = if (large_base_size) { 22 } else { 12 }))
  if (log2_transform) {
    g <- g + ylab('log2 count')
  }
  if (grey) {
    g <- g + theme(panel.background = element_rect(fill = '#f2f2f2',
                                                   colour = '#f2f2f2',
                                                   size = 0.5, linetype = 'solid'))
  } else {
    g <- g + theme_minimal()
  }
  
  list(plot              = g,
       stats             = cbind(stats, inds_char = 'NULL'),
       pseudotime_bounds = pseudotime_bounds)
}


