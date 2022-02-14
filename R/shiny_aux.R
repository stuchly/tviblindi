#### Shiny interface for tviblindi: auxiliary functions

## Function: remove duplicit rows from data frame
.row_unique <- function(df) {
    hashes <- apply(df, 1, function(x) paste(x, collapse = ' '))
    dupes  <- which(duplicated(hashes))
    if (length(dupes) > 0) return(df[-dupes, ])
    df
}

## Function: list of numeric vectors -> list of strings
.label <- function(l, collapse = ' ') {
    lapply(l, function(x) {
        out <- paste(x, collapse = collapse)
        if (out == '') return('nil')
        out
    })
}

## Function: find slot in list containing all elements of given vector
.find <- function(what, where) which(sapply(where, function(x) all(what %in% x)))

## Function: hclust reference -> hclust leaf node (recursive)
.ref_to_vals <- function(ref, where) {

    vals <- vector(mode = 'integer')
    vec  <- where[ref, ]
    vals <- c(vals, -vec[vec < 0])
    for (i in vec[vec > 0]) vals <- c(vals, .ref_to_vals(i, where))
    vals
}

trajectories_dendrogram <- function(precomputed_dendrogram         = NULL,
                                    precomputed_dendrogram_labels  = NULL,
                                    zoom_idcs                      = NULL,
                                    leaves_to_highlight.A          = NULL,
                                    leaves_to_highlight.B          = NULL,
                                    leaves_to_highlight.zoom       = NULL,
                                    pers                           = NULL,
                                    repre.reduced                  = NULL,
                                    perc                           = NULL,
                                    out.data                       = NULL,
                                    out.classif                    = NULL,
                                    out.labels                     = NULL,
                                    out.dendrogram                 = NULL,
                                    make_png_dendrogram            = FALSE) {
    ## Create ggplot2 trajectories dendrogram & (optionally) output an hclust object
    if (is.null(precomputed_dendrogram)) {
        if (perc == 100) return(FALSE)
        R              <- repre.reduced; rm(repre.reduced)
        Ru             <- c(unique(unlist(R)))            # unique simplices

        threshold      <- ceiling(length(R) / 100 * perc) # lower bound of num trajectories per leaf
        H <- lapply(Ru, function(u) {                     # idcs of simplices corresponding to deaths of relevant homology classes
            idcs <- which(pers$inds$death == u)
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
        M                        <- matrix(M, nc = 2, byrow = TRUE)

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
        data              <- dendro_data(d, type = 'rectangle')
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
        data              <- dendro_data(d, type = 'rectangle')
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

## Function: resolve binary tree branching (recursive)
.recurse <- function(pres, idcs, n, father, son, brother, th = 1) {

    if (length(idcs) > th & n > 0) {
        left <- which(unlist(lapply(idcs, FUN = function(i) pres$Ru[n] %in% pres$R[[i]])))

        while (length(left) == 0 || length(left) == length(idcs)) { # only divide if purity is <1
            n    <- n - 1
            left <- which(unlist(lapply(idcs, FUN = function(i) pres$Ru[n] %in% pres$R[[i]])))

            if (n == 0) {
                pres$leafcount              <- pres$leafcount + 1
                pres$leaves[pres$leafcount] <- son
                col                         <- list(n,
                                                    length(idcs),
                                                    idcs,
                                                    father,
                                                    NULL,
                                                    NULL,
                                                    brother)
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

        col  <- list(n,
                     length(idcs),
                     NULL,
                     father,
                     sonl,
                     sonr,
                     brother)

        assign(son, col, envir = pres)

        .recurse(pres, idcs[left],  n - 1, son, sonl, brother, th)
        .recurse(pres, idcs[-left], n - 1, son, sonr, brother, th)
    } else { # length(idxs) <= th | n <= 0
        pres$leafcount              <- pres$leafcount + 1
        pres$leaves[pres$leafcount] <- son
        col                         <- list(n,
                                            length(idcs),
                                            idcs,
                                            father,
                                            NULL,
                                            NULL,
                                            brother)
        assign(son, col, envir = pres)
    }
}

.draw_placeholder <- function(picture = 'tree') {
    if (picture == 'tree') {
        j   <- jpeg::readJPEG(system.file('tree.jpg', package = 'tviblindi'), native = TRUE)
    }
    if (picture == 'petal') {
        j   <- jpeg::readJPEG(system.file('petal.jpg', package = 'tviblindi'), native = TRUE)
    }
    plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
    rasterImage(j, 0, 0, 1, 1)
}

.add_path_info_to_fcs <- function(ff, pseudotime, all_walks, trajectories_of_interest, event_sel = NULL, id = NULL) {
    pp            <- unique(select_paths_points(all_walks, trajectories_of_interest))
    out           <- matrix(-1, nrow = nrow(ff@exprs), ncol = 2)
    colnames(out) <- c(paste(id,'which_event', sep = '_'), paste(id, 'local_pseudotime', sep = '_'))
    if (is.null(event_sel)) event_sel <- 1:nrow(ff@exprs)

    out[event_sel[pp], 1] <- 1000
    out[event_sel[pp], 2] <- as.numeric(as.factor(pseudotime$res[pp]))

    make_valid_fcs(cbind(ff@exprs, out), desc1 = as.character(ff@parameters@data$desc))
}

## Function: add a column (~ channel) to an FCS file
fcs.add_col <- function(ff, new_col, colname = 'label') {
    require(flowCore)

    if (class(ff) != 'flowFrame') stop('Parameter ff is not a valid flowFrame object')

    efcs <- ff@exprs
    N    <- nrow(efcs)
    len  <- length(new_col)

    if (N != len)  stop(paste0('Number of rows of expression matrix is ', N, ', whereas length of new column is ', len, '.'))

    params <- ff@parameters
    pd     <- pData(params)
    cols   <- as.vector(pd$name)
    idcs   <- match(cols, pd$name)

    if (any(is.na(idcs))) stop('Invalid column specifier')

    channel_number     <- ncol(ff) + 1
    channel_id         <- paste0('$P', channel_number)
    channel_name       <- colname
    channel_range      <- max(new_col) + 1
    channel_min        <- min(0, min(new_col) - 1)
    plist              <- matrix(c(channel_name, channel_name, channel_range,
                                   channel_min, channel_range - 1))
    rownames(plist)    <- c('name', 'desc', 'range', 'minRange', 'maxRange')
    colnames(plist)    <- c(channel_id)
    pd                 <- rbind(pd, t(plist))
    pData(params)      <- pd
    channel_names      <- colnames(efcs)
    efcs.mod           <- cbind(efcs, new_col)
    colnames(efcs.mod) <- c(channel_names, colname)

    ff.mod             <- flowFrame(efcs.mod, params, description = keyword(ff))

    keyval                                      <- list()
    keyval[[paste0('$P', channel_number, 'B')]] <- '32'
    keyval[[paste0('$P', channel_number, 'R')]] <- toString(channel_range)
    keyval[[paste0('$P', channel_number, 'E')]] <- '0,0'
    keyval[[paste0('$P', channel_number, 'N')]] <- channel_name
    keyval[[paste0('$P', channel_number, 'S')]] <- channel_name
    keyword(ff.mod)                             <- keyval

    flowCoreP_Rmax <- paste0('flowCore_$P', channel_number, 'Rmax')
    flowCoreP_Rmin <- paste0('flowCore_$P', channel_number, 'Rmin')

    keyword(ff.mod)[flowCoreP_Rmax] <- max(20000, keyword(ff.mod)$`flowCore_$P1Rmax`)
    keyword(ff.mod)[flowCoreP_Rmin] <- 0

    ff.mod
}

.plot_trajectories <- function(X,
                               walks,
                               walk_idcs.A,
                               walk_idcs.B,
                               flip_colours,
                               pseudotime_highlight_bounds,
                               pseudotime,
                               highlight_in_background,
                               selected_trajectory_points=NULL,
                               pointsize = .05,
                               ...) {
    ## Plot trajectories over a 2D layout
    col1 <- c(34, 87, 201, 255)
    col2 <- c(194, 45, 55, 255)
    plot(scattermore(X, rgba = c(200, 200, 200, 150), cex = pointsize))

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

            plot(scattermore(pts, rgba = c(0, 153, 31, 255), xlim = c(0, 1), ylim = c(0, 1), cex = pointsize + 1), add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
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
        plot(scattermore(pts1, rgba = col1, xlim = c(0, 1), ylim = c(0, 1)), add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
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
        plot(scattermore(pts2, rgba = col2, xlim = c(0, 1), ylim = c(0, 1)), add = TRUE)
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
            plot(scattermore(pts, rgba = c(0, 153, 31, 255), xlim = c(0, 1), ylim = c(0, 1), cex = pointsize + 1), add = TRUE, xlim = c(0, 1), ylim = c(0, 1))
        }
    }
}

.plot_trajectories_full <- function(X,
                                    walks,
                                    walk_idcs.A,
                                    walk_idcs.B,
                                    flip_colours,
                                    pseudotime_highlight_bounds,
                                    pseudotime,
                                    highlight_in_background,
                                    pointsize = .08,
                                    ...) {
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
        lines(pts, lty = 1, col = alpha('darkgrey', 0.02), ...)
        points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = alpha('darkgrey', 0.02))
    }
    if (length(sel1) > 0) {
        for (i in sel1) {
            j   <- j + 1
            pts <- X[i, ]
            lines(pts, lty = 1, col = alpha(cols[1], 0.3), ...)
            points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = alpha(cols[1], 0.3))
        }
    }
    if (length(sel2) > 0) {
        for (i in sel2) {
            j   <- j + 1
            pts <- X[i, ]
            lines(pts, lty = 1, col = alpha(cols[2], 0.3), ...)
            points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = alpha(cols[2], 0.3))
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
.plot_tracked_markers <- function(walks,
                                  walk_idcs,
                                  tv,
                                  pseudotime,
                                  markers,
                                  breaks   = NULL,
                                  n.part   = 10,
                                  exp.part = 1,
                                  large_base_size = FALSE,
                                  grey = FALSE) {
    require(ggplot2)

    if (length(markers) == 1) {
        return(do.call(.plot_tracked_markers_single,
                       as.list(environment())))
    } else {
        return(do.call(.plot_tracked_markers_multiple,
                       as.list(environment())))
    }
}

.plot_tracked_markers_single <- function(walks,
                                         walk_idcs,
                                         tv,
                                         pseudotime,
                                         markers,
                                         breaks      = NULL,
                                         n.part      = 10,
                                         exp.part    = 1,
                                         show_legend = FALSE,
                                         large_base_size = FALSE,
                                         grey            = FALSE) {
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
        b <- as.numeric(quantile(pp, p))
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
                                        # for each point on walk, which segment does it fall into?

    pseudotime_bounds <- sort(c(unlist(progress)[which(!duplicated(categs))], 1))
    coords <- tv$data[, markers]

    stats  <- lapply(1:N, function(i) {
        inds <- categs == i # pick points on paths by pseudotime increment
        if (!any(inds)) return(NULL)

        which.walks <- unlist(lapply(1:length(walks), function(j) rep(j, length(walks[[j]]))))[inds]
        pts         <- unlist(walks)[unlist(inds)]
        means       <- lapply(sort(unique(which.walks)), function(j) mean(coords[pts[which.walks == j]]))

        inds_char  <- lapply(sort(unique(which.walks)), function(j) paste(pts[which.walks == j], collapse =","))

        ## v           <- cbind(rep(i, length(means)), sort(unique(which.walks)), unlist(means), unlist(inds_char))
        ## if (is.null(dim(v))) { names(v) <- c('segment', 'walk', 'expression', 'inds_char') } else { colnames(v) <- c('segment', 'walk', 'expression','inds_char') }

        v           <- cbind(rep(i, length(means)), sort(unique(which.walks)), unlist(means))
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

.plot_tracked_markers_multiple <- function(walks,
                                           walk_idcs,
                                           tv,
                                           pseudotime,
                                           markers,
                                           breaks      = NULL,
                                           n.part      = 10,
                                           exp.part    = 1,
                                           show_legend = TRUE,
                                           large_base_size = FALSE,
                                           grey = FALSE) {
    pseudotime <- pseudotime$res
    pseudotime <- pseudotime/max(pseudotime)

    walks    <- walks[walk_idcs]
    progress <- lapply(walks, function(pts) pseudotime[pts])
    pp       <- unique(unlist(progress))

    b <- NULL
    N <- NULL

    if (is.null(breaks)) {
        p <- (seq(0, 1, 1 / (n.part))) ^ exp.part
        b <- as.numeric(quantile(pp, p))
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
    list(plot              = g,
         stats             = cbind(stats,inds_char="NULL"),
         pseudotime_bounds = pseudotime_bounds)
}

.plot_tracked_populations <- function(walks,
                                      walk_idcs,
                                      tv,
                                      labels_name,
                                      pseudotime,
                                      populations,
                                      breaks = NULL,
                                      n.part = 20,
                                      exp.part = 1,
                                      large_base_size = FALSE,
                                      log2_transform,
                                      grey = FALSE) {
    pseudotime <- pseudotime$res
    pseudotime <- pseudotime / max(pseudotime)

    walks    <- walks[walk_idcs]
    progress <- lapply(walks, function(pts) pseudotime[pts])
    pp       <- unique(unlist(progress))

    b <- NULL
    N <- NULL

    if (is.null(breaks)) {
        p <- (seq(0, 1, 1 / (n.part))) ^ exp.part
        b <- as.numeric(quantile(pp, p))
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

.update_walks_by_termini <- function(tv,
                                     pathmodel_name,
                                     pseudotime,
                                     marked_termini,
                                     termini_per_path,
                                     death_birth_ratio,
                                     death_on_x_axis) {
    ## Identify chosen walks
    idcs <- which(termini_per_path %in% marked_termini)

    ## Remove other walks
    walks.selected <- lapply(idcs, function(idx) select_paths_points(tv$walks[[pathmodel_name]], idx))
    lens           <- sapply(walks.selected, length)
    walks.selected <- list(v      = unlist(walks.selected),
                           starts = c(1, 1 + cumsum(lens[-length(lens)])))

    ## Put terminal node with highest pseudotime at end of each walk
    if (length(marked_termini) > 1) {
        max_t             <- marked_termini[which.max(pseudotime$res[marked_termini])][1]
        ends              <- c(walks.selected$starts[-1] - 1, length(walks.selected$v))
        walks.selected$v[ends] <- max_t
    }

    ## Cluster and triangulate walks
    withProgress(message = 'Contracting trajectories', expr = {
        walks_clusters <- remove_cycles(contract_walks(walks.selected, tv$clusters), verbose = FALSE)

        sel          <- 1:length(walks_clusters$starts)
        s2           <- which(unlist(lapply(tv$filtration$cmplx, function(x) length(x) == 2))) # 2-simplex idcs
        cmplx        <- tv$filtration$cmplx[s2]                                        # 2-simplices
        cmplx_hashed <- hash_cmplx(cmplx)                                              # hash the simplices list for faster look-up
        edges        <- do.call(rbind, lapply(cmplx, function(sim, XX) return(c(sim,
                                                                                dist(tv$codes[sim,]))),
                                              XX = tv$codes))
        colnames(edges) <-c('from', 'to', 'weight')

        graph           <- igraph::graph_from_edgelist(edges[, 1:2])
        E(graph)$weight <- edges[, 3]
        triangulation   <- list()[1:length(sel)]

        j <- 0
    })

    withProgress(message = 'Triangulating paths', expr = {
        for (i in sel) {
            j <- j + 1
            triangulation[[j]] <- s2[onepath_triangulation(select_paths_points(walks_clusters, i),
                                                           tv$codes, cmplx_hashed, graph)]
        }
    })

    N <- length(idcs)
    repre      <- list()[1:N]
    repre[[1]] <- integer(0)
    j          <- 1

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


    ##METHOD CHANGED
    ss<-which(unlist(lapply(repre[-1],function(x) any(x<0))))+1
    if (length(ss)>0){
        print(length(ss))
        withProgress(message = 'Recalculating Inf value representations - just this once', value = tick, expr = {
            for (idx in ss) {
                print(idx)
                cycle      <- path_sum(triangulation[[1]], triangulation[[idx]])
                this_repre <- get_rep_straight(cycle, tv$reduced_boundary, tv$boundary, update = TRUE)
                repre[[idx]] <- this_repre
                incProgress(tick)
            }
        })

    }

    walks <- lapply(1:N, function(idx) { select_paths_points(walks.selected, idx) })


    p <- .compute_persistence(tv, repre, death_birth_ratio, death_on_x_axis)

    return(list(random_walks = walks,
                repre        = repre,
                pers         = p$pers,
                pers_diag    = p$pd,
                triangulation=triangulation,
                repre=repre))
}

##METHOD CHANGED
.update_walks_by_dendrogram <-  function(tv,
                                         pathmodel_name,
                                         pseudotime,
                                         marked_termini,
                                         idcs,
                                         walks_sel,
                                         triangulation,
                                         termini_per_path,
                                         death_birth_ratio,
                                         death_on_x_axis) {
    ## Identify chosen walks
    triangulation<-triangulation[idcs]

    idcs <- walks_sel[idcs]

    ## Remove other walks
    walks.selected <- lapply(idcs, function(idx) select_paths_points(tv$walks[[pathmodel_name]], idx))
    lens           <- sapply(walks.selected, length)
    walks.selected <- list(v      = unlist(walks.selected),
                           starts = c(1, 1 + cumsum(lens[-length(lens)])))

    ## Put terminal node with highest pseudotime at end of each walk
    if (length(marked_termini) > 1) {
        max_t             <- marked_termini[which.max(pseudotime$res[marked_termini])][1]
        ends              <- c(walks.selected$starts[-1] - 1, length(walks.selected$v))
        walks.selected$v[ends] <- max_t
    }

    ## ## Cluster and triangulate walks
    ## withProgress(message = 'Contracting trajectories', expr = {
    ##     walks_clusters <- remove_cycles(contract_walks(walks.selected, tv$clusters), verbose = FALSE)

    ##     sel          <- 1:length(walks_clusters$starts)
    ##     s2           <- which(unlist(lapply(tv$filtration$cmplx, function(x) length(x) == 2))) # 2-simplex idcs
    ##     cmplx        <- tv$filtration$cmplx[s2]                                        # 2-simplices
    ##     cmplx_hashed <- hash_cmplx(cmplx)                                              # hash the simplices list for faster look-up
    ##     edges        <- do.call(rbind, lapply(cmplx, function(sim, XX) return(c(sim,
    ##                                                                             dist(tv$codes[sim,]))),
    ##                                           XX = tv$codes))
    ##     colnames(edges) <-c('from', 'to', 'weight')

    ##     graph           <- igraph::graph_from_edgelist(edges[, 1:2])
    ##     E(graph)$weight <- edges[, 3]
    ##     triangulation   <- list()[1:length(sel)]

    ##     j <- 0
    ## })

    ## withProgress(message = 'Triangulating paths', expr = {
    ##     for (i in sel) {
    ##         j <- j + 1
    ##         triangulation[[j]] <- s2[onepath_triangulation(select_paths_points(walks_clusters, i),
    ##                                                        tv$codes, cmplx_hashed, graph)]
    ##     }
    ## })

    N <- length(idcs)
    repre      <- list()[1:N]
    repre[[1]] <- integer(0)
    j          <- 1

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
    walks <- lapply(1:N, function(idx) { select_paths_points(walks.selected, idx) })

    p <- .compute_persistence(tv, repre, death_birth_ratio, death_on_x_axis)

    return(list(random_walks = walks,
                repre        = repre,
                pers         = p$pers,
                pers_diag    = p$pd,
                walks_sel=idcs,
                triangulation=triangulation))
}


.compute_persistence <- function(
                                 tv,
                                 repre,
                                 death_birth_ratio,
                                 death_on_x_axis
                                 ) {
    pers <- pers_diagram(dBr = tv$reduced_boundary, repre = repre, plot = FALSE)
    ##METHOD CHANGED
    if(any(is.infinite(pers$vals$death))){

        ##print("Yes!!")
        ss<-which(is.infinite(pers$vals$death))
        ##print(-max(pers$vals$death[-ss]))
        m_ypos<-max((pers$vals$death[-ss]-pers$vals$birth[-ss])/2)
        m_ssb<-pers$vals$birth[ss]
        #print(m_ypos)
        ## print(str(tv1$reduced_boundary))
        ## print(str(tv1$boundary))
        pers$vals$death[ss]<-(2.1*m_ypos+m_ssb)
        message(length(ss), " essential homology class(es) detected.")
        message("position: x=(",paste((pers$vals$death[ss]+pers$vals$birth[ss])/2,collapse=","),")\n assigned death=(",paste(pers$vals$death[ss],collapse=","),"\n assigned (death-birth)/2 =(",paste((pers$vals$death[ss]-pers$vals$birth[ss])/2,collapse=","))
        ##ss1<-which.max(pers$vals$death[-ss])
        ##print((pers$vals$birth[-ss][ss1]+pers$vals$death[-ss][ss1])/2)
        ##print((pers$vals$death[ss]+pers$vals$birth[ss])/2)
    }
    pd <- .persistence_diagram(pers, death_birth_ratio, death_on_x_axis)
    list(
        pers = pers,
        pd = pd
    )
}

.persistence_diagram <- function(
                                 pers,
                                 death_birth_ratio,
                                 death_on_x_axis
                                 ) {
    pd   <- data.frame(Dimension    = pers$vals$dim,
                       Birth        = pers$vals$birth,
                       Death        = pers$vals$death,
                       BirthSimplex = pers$inds$birth,
                       DeathSimplex = pers$inds$death,
                       xplot        = if (death_on_x_axis) { pers$vals$death } else { (pers$vals$birth + pers$vals$death) / 2 },
                       yplot        = if (death_birth_ratio) { pers$vals$death / pers$vals$birth } else { (pers$vals$death - pers$vals$birth) / 2 })
    pd[pd$Dimension > 0, ]
}