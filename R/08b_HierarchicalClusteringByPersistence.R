.row_unique <- function(
  ## Remove duplicit rows from data frame
  df
) {
  hashes <- apply(df, 1, function(x) paste(x, collapse = ' '))
  dupes  <- which(duplicated(hashes))
  if (length(dupes) > 0) return(df[-dupes, ])
  df
}


.label <- function(
  ## List of numeric vectors -> list of strings
  l,
  collapse = ' '
) {
  lapply(l, function(x) {
    out <- paste(x, collapse = collapse)
    if (out == '') return('nil')
    out
  })
}


.find <- function(
  ## Find slot in list containing all elements of given vector
  what,
  where
) {
  which(
    sapply(
      where,
      function(x)
        all(what %in% x)
    )
  )
}

.ref_to_vals <- function(
  ## hclust reference -> hclust leaf node (recursive)
  ref,
  where
) {
  vals <- vector(mode = 'integer')
  vec  <- where[ref, ]
  vals <- c(vals, -vec[vec < 0])
  for (i in vec[vec > 0]) vals <- c(vals, .ref_to_vals(i, where))
  vals
}

trajectories_dendrogram <- function(
  precomputed_dendrogram        = NULL,
  precomputed_dendrogram_labels = NULL,
  zoom_idcs                     = NULL,
  leaves_to_highlight.A         = NULL,
  leaves_to_highlight.B         = NULL,
  leaves_to_highlight.zoom      = NULL,
  pers                          = NULL,
  repre.reduced                 = NULL,
  perc                          = NULL,
  out.data                      = NULL,
  out.classif                   = NULL,
  out.labels                    = NULL,
  out.dendrogram                = NULL,
  make_png_dendrogram           = FALSE
) {
  ## Prevent 'no visible binding for global variable' build warning
  x <- y <- xend <- yend <- label <- NULL
  
  ## Create ggplot2 trajectories dendrogram & (optionally) output an hclust object
  if (is.null(precomputed_dendrogram)) {
    if (perc == 100) return(FALSE)
    R              <- repre.reduced; rm(repre.reduced)
    Ru             <- c(unique(unlist(R)))            # unique simplices
    
    threshold      <- ceiling(length(R) / 100 * perc) # lower bound of num trajectories per leaf
    H <- lapply(Ru, function(u) {                     # idcs of simplices corresponding to deaths of relevant homology classes
      idcs <- which(pers$idcs$death == u)
      if (length(idcs) == 0) return()
      idcs
    })
    H              <- unlist(H)
    H.vals         <- sort(pers$vals$death[H])        # filtration values associated with deaths of relevant homology classes
    Ru             <- Ru[order(H)]                    # sort simplex idcs by corresponding filtration values
    Rl             <- .label(R)                       # string hashes for reduced trajectory representations
    R.unique_idcs  <- sapply(unique(Rl), function(r) which(Rl == r)[1])
    # idcs of unique reduced trajectory representations
    R.unique       <- R[R.unique_idcs]                # only keep those unique representations
    Rl.unique      <- Rl[R.unique_idcs]
    
    ## Build binary tree
    M              <- do.call(rbind, lapply(R.unique, function(r) Ru %in% r))
    N              <- length(Ru)
    K              <- length(Rl.unique)
    
    pres           <- new.env(hash = TRUE)
    pres$R         <- R
    pres$Ru        <- Ru
    pres$ind       <- 0
    pres$leafcount <- 0
    pres$leaves    <- list()
    pres$mcount    <- 0
    pres$merges    <- list()
    pres$hcount    <- 0
    pres$heights   <- list()
    
    ## Establish binary tree branching recursively
    .recurse(pres, 1:length(pres$R), length(pres$Ru), 'NULL', 'ROOT', 'NULL', th = threshold)
    
    leaves  <- pres$leaves
    h       <- unlist(pres$heights)
    m       <- do.call(rbind, pres$merges)[order(h), , drop = FALSE]
    h       <- sort(h)
    
    if (is.null(dim(m))) return(FALSE)
    
    ## Create hclust merging pattern
    l                        <- lapply(leaves, function(leaf) which(apply(m, 1, function(x) leaf %in% x))); names(l) <- leaves
    mm                       <- as.vector(t(m))
    M                        <- rep(0, length(mm))
    M[which(mm %in% leaves)] <- -(1:length(leaves))
    M                        <- matrix(M, ncol = 2, byrow = TRUE)
    
    for (i in 1:nrow(M)) {
      for (j in 1:2) {
        if (M[i, j] == 0) {
          node      <- m[i, j]
          child     <- get(node, pres)[[5]]
          M[i, j]   <- l[[child]]
          l[[node]] <- i
        }
      }
    }
    
    merge <- M; rm(M)
    
    ## Resolve leaf ordering for a valid dendrogram
    o <- as.list(1:length(leaves))
    for (i in 1:nrow(merge)) {
      
      p <- merge[i, ]
      
      if (sum(p < 0) == 2) {
        
        idx1         <- .find(abs(p[1]), o)
        idx2         <- .find(abs(p[2]), o)
        tmp          <- o[[idx2]]
        o[[idx2]]    <- integer(0)
        o[[idx1]]    <- c(o[[idx1]], tmp)
        
      } else if (sum(p < 0) == 1) {
        
        which.ref    <- which(p > 0)
        ref          <- p[which.ref]
        val          <- -p[3 - which.ref]
        idx.ref      <- .find(.ref_to_vals(p[which.ref], merge), o)
        idx.val      <- .find(val, o)
        tmp          <- o[[idx.val]]
        o[[idx.val]] <- integer(0)
        o[[idx.ref]] <- c(o[[idx.ref]], tmp)
        
      } else if (sum(p < 0) == 0) {
        
        idx1         <- .find(.ref_to_vals(p[1], merge), o)
        idx2         <- .find(.ref_to_vals(p[2], merge), o)
        tmp          <- o[[idx2]]
        o[[idx2]]    <- integer(0)
        o[[idx1]]    <- c(o[[idx1]], tmp)
      }
    }
    
    ## Create a valid hclust object
    o           <- unlist(o)
    tree        <- list()
    tree$merge  <- merge
    tree$order  <- o
    tree$height <- H.vals[h]
    nodes       <- as.vector(t(m))
    labs        <- nodes[nodes %in% leaves]
    tree$labels <- labs
    class(tree) <- 'hclust'
    
    d                 <- as.dendrogram(tree)
    data              <- ggdendro::dendro_data(d, type = 'rectangle')
    labs              <- as.character(data$labels$label)
    tree_data         <- data
    leaves            <- labs
    labs_per_leaf     <- as.character(sapply(as.character(data$labels$label), function(leaf) get(leaf, pres)[[2]]))
    
    if (!is.null(out.data)) eval.parent(substitute(out.data <- tree_data))
    data$labels$label <- labs_per_leaf
    
    if (!is.null(out.labels)) eval.parent(substitute(out.labels <- list(labs,
                                                                        labs_per_leaf)))
    if (!is.null(out.dendrogram)) eval.parent(substitute(out.dendrogram <- d))
    
    cl          <- lapply(labs, function(leaf) get(leaf, pres)[[3]]) # walk indices per leaf
    names(cl)   <- labs
    if (!is.null(zoom_idcs)) {
      cl <- cl[zoom_idcs[1]:zoom_idcs[2]]
    }
    if (!is.null(out.classif)) eval.parent(substitute(out.classif <- cl))
    
  } else { # !is.null(precomputed_dendrogram)
    d             <- precomputed_dendrogram
    labs          <- precomputed_dendrogram_labels[[1]]
    labs_per_leaf <- precomputed_dendrogram_labels[[2]]
    data              <- ggdendro::dendro_data(d, type = 'rectangle')
    tree_data         <- data
    leaves            <- labs
    if (!is.null(out.data)) eval.parent(substitute(out.data <- tree_data))
    data$labels$label <- labs_per_leaf
  }
  
  if (!is.null(zoom_idcs)) {
    X            <- c(data$labels[zoom_idcs, 1])
    x_min        <- min(X)
    x_max        <- max(X)
    # remove segments not in zoom selection
    which_remove <- which(apply(data$segments, 1, function(v)
      (min(v[c(1, 3)]) < x_min || min(v[c(1, 3)]) > x_max) &&
        (max(v[c(1, 3)]) > x_max || max(v[c(1, 3)]) < x_min)
    ))
    # cut under-/overflowing segments short
    if (length(which_remove) > 0) { data$segments <- data$segments[-which_remove, ] }
    overflow  <- which(data$segments[, 1] > x_max)
    underflow <- which(data$segments[, 1] < x_min)
    data$segments[overflow, 1]  <- x_max
    data$segments[underflow, 1] <- x_min
    overflow  <- which(data$segments[, 3] > x_max)
    underflow <- which(data$segments[, 3] < x_min)
    data$segments[overflow, 3]  <- x_max
    data$segments[underflow, 3] <- x_min
    # remove labels not in zoom selection
    which_remove <- which(data$labels[, 1] > x_max | data$labels[, 1] < x_min)
    if (length(which_remove) > 0) { data$labels <- data$labels[-which_remove, ] }
  }
  
  p.png <- NULL
  
  if (make_png_dendrogram) {
    p.png <- ggplot(data$segments) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), lineend = 'round', linejoin = 'round',
                                                  size = if (!is.null(zoom_idcs)) { 2 } else { 1 }) +
      geom_text(data = data$labels, aes(x, y - 0.01, label = label), hjust = 1, angle = 0, size = if (!is.null(zoom_idcs)) { 6 } else { 4.1 }) +
      cowplot::theme_nothing() +
      theme(plot.margin = unit(c(-.2, -0.05, -.2, 0), 'cm')) +
      coord_flip()
  }
  
  p.regular <- ggplot(data$segments) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend), lineend = 'round', linejoin = 'round',
                                                    size = if (!is.null(zoom_idcs)) { .8 } else { .5 }) +
    geom_text(data = data$labels, aes(x, y, label = label), hjust = 1, angle = 0, size = if (!is.null(zoom_idcs)) { 5.2 } else { 3.4 }) +
    cowplot::theme_nothing() +
    theme(plot.margin = unit(c(-.2, 0, -.2, 0), 'cm')) +
    coord_flip()
  
  if (!is.null(leaves_to_highlight.zoom) || (!is.null(leaves_to_highlight.A) || !is.null(leaves_to_highlight.B))) {
    branches           <- data.frame(tree_data$segments[tree_data$segments$yend == 0, ], label = tree_data$labels$label)
    colnames(branches) <- c('xmin', 'ymin', 'xmax', 'ymax', 'label')
    ymax <- max(tree_data$segments[, c(2, 4)])
  }
  
  if (!is.null(leaves_to_highlight.zoom)) {
    lowerbound <- min(branches$xmin[leaves_to_highlight.zoom])
    upperbound <- max(branches$xmin[leaves_to_highlight.zoom])
    rect <- geom_rect(xmin = lowerbound - 0.2, xmax = upperbound + 0.2, ymin = ymax *.96, ymax = ymax, fill = '#9cf0d9', alpha = .01)
    p.regular <- p.regular + rect
    if (make_png_dendrogram) {
      p.png <- p.png + rect
    }
  }
  
  if (!is.null(leaves_to_highlight.A) || !is.null(leaves_to_highlight.B)) {
    
    divide_leaves_by_subtrees <- function(leaf_idcs) {
      if (length(leaf_idcs) == 1) return(list(leaf_idcs))
      subtrees <- vector(mode = 'list')
      n_leaves <- 0
      tmp <- c(leaf_idcs[1])
      idx <- 2
      while (idx <= length(leaf_idcs)) {
        this <- leaf_idcs[idx]
        if (this == tmp[length(tmp)] + 1) {
          tmp <- c(tmp, this)
          if (idx == length(leaf_idcs)) {
            subtrees[[length(subtrees) + 1]] <- tmp
          }
        } else {
          subtrees[[length(subtrees) + 1]] <- tmp
          tmp <- c(this)
        }
        idx <- idx + 1
      }
      return(subtrees)
    }
    
    if (!is.null(leaves_to_highlight.A)) {
      idcs     <- which(branches$label %in% leaves_to_highlight.A)
      X        <- branches$xmin[idcs]
      subtrees <- divide_leaves_by_subtrees(X)
      for (X in subtrees) {
        lowerbound <- min(X)
        upperbound <- max(X)
        rect <- geom_rect(xmin = lowerbound - 0.25, xmax = upperbound + 0.25, ymin = 0, ymax = ymax * .95, fill = '#c2daff',
                          alpha = if (!is.null(zoom_idcs)) { .03 } else { .01 })
        p.regular <- p.regular + rect
        if (make_png_dendrogram) {
          p.png <- p.png + rect
        }
      }
    }
    if (!is.null(leaves_to_highlight.B)) {
      idcs     <- which(branches$label %in% leaves_to_highlight.B)
      X        <- branches$xmin[idcs]
      subtrees <- divide_leaves_by_subtrees(X)
      for (X in subtrees) {
        lowerbound <- min(X)
        upperbound <- max(X)
        rect <- geom_rect(xmin = lowerbound - 0.2, xmax = upperbound + 0.2, ymin = 0, ymax = ymax * .95, fill = 'pink',
                          alpha = if (!is.null(zoom_idcs)) { .05 } else { .01 })
        p.regular <- p.regular + rect
        if (make_png_dendrogram) {
          p.png <- p.png + rect
        }
      }
    }
  }
  
  list(
    regular = p.regular,
    png     = p.png
  )
}

.recurse <- function(
  ## Resolve binary tree branching (recursive)
  pres,
  idcs,
  n,
  father,
  son,
  brother,
  th = 1
) {
  if (length(idcs) > th & n > 0) {
    left <- which(unlist(lapply(idcs, FUN = function(i) pres$Ru[n] %in% pres$R[[i]])))
    
    while (length(left) == 0 || length(left) == length(idcs)) { # only divide if purity is <1
      n    <- n - 1
      left <- which(unlist(lapply(idcs, FUN = function(i) pres$Ru[n] %in% pres$R[[i]])))
      
      if (n == 0) {
        pres$leafcount              <- pres$leafcount + 1
        pres$leaves[pres$leafcount] <- son
        col <- list(
          n,
          length(idcs),
          idcs,
          father,
          NULL,
          NULL,
          brother
        )
        assign(son, col, envir = pres)
        return(-1)
      }
    }
    sonl                        <- paste(son,  n,  sep = '_')
    sonr                        <- paste(son, -n, sep = '_')
    
    pres$mcount                 <- pres$mcount + 1
    pres$merges[[pres$mcount]]  <- c(sonl, sonr)
    pres$hcount                 <- pres$hcount + 1
    pres$heights[[pres$hcount]] <- n
    
    col <- list(
      n,
      length(idcs),
      NULL,
      father,
      sonl,
      sonr,
      brother
    )
    
    assign(son, col, envir = pres)
    
    .recurse(pres, idcs[left],  n - 1, son, sonl, brother, th)
    .recurse(pres, idcs[-left], n - 1, son, sonr, brother, th)
  } else { # length(idxs) <= th | n <= 0
    pres$leafcount              <- pres$leafcount + 1
    pres$leaves[pres$leafcount] <- son
    col <- list(
      n,
      length(idcs),
      idcs,
      father,
      NULL,
      NULL,
      brother
    )
    assign(son, col, envir = pres)
  }
}

