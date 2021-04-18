#' Constructor of a tviblindi-analysis object
#'
#' Constructs an object of class \code{tviblindi} for trajectory-inference analysis.
#' This object allows for semi-automated detection of developmental trajectories in expression data acquired via high-dimensional cytometry, scRNA-seq or CITE-seq ADT data.
#' 
#' After calling this constructor, use the methods \code{SetOrigin}, \code{ConstructkNNG}, \code{Denoise}, \code{AddLayout}, \code{Cluster}, \code{Filter}, \code{ComputePseudotime}, \code{SimulateRandomWalks} and \code{Interactive} to conduct you trajetory-inference analysis.
#' 
#' Alternatively, use the function \code{InferTrajectories} instead of using this constructor and these methods to run the whole \code{tviblindi} pipeline.
#'
#' @param data numeric matrix: pre-processed protein expression or read count data matrix. Column names shall correspond to name of individual markers (detected on a transctipt or protein level)
#' @param labels string or factor vector of the same length as row count of \code{data}. Assigns a label to each measured event, mapping the event to a known cell population. Unassigned cells should be labelled as '\code{ungated}' or '\code{\*ungated\*}'
#' If \code{labels} is left empty, all events are assigned the '\code{ungated}' label
#' @param fcs_path optional string: path to FCS file from which cytometry expression \code{data} was extracted. The use of \code{fcs_path} is encouraged if a subset of expression data from a known FCS file is used. If \code{fcs_path} is available, an enhanced version of the FCS file (with information about detected developmental pathways) can be generated for downstream analysis in FlowJo. Default value is \code{NULL}
#' @param fcs_subset_idcs integer vector: row numbers in the expression matrix of an associated FCS file that correspond to rows of \code{data}.
#' (Default value is \code{NULL}, denoting that all events from the associated FCS file, if any such file exists, were taken along for the trajectory-inference analysis)
#' @param analysis_name optional string: name of the trajectory-inference analysis. If not specified, a name is created by concatenation of user-name, date, time and \code{tviblindi} version
#'
#' @return \code{tviblindi} returns an object of class \code{tviblindi}
#'
#' @export
tviblindi <- function(
  data,
  labels,
  fcs_path = NULL,
  fcs_subset_idcs = NULL,
  analysis_name = NULL
) {
  
  if (is.null(analysis_name)) {
    username <- Sys.info()['login']
    if (is.null(username))
      username <- Sys.info()['user']
    if (is.null(username))
      username <- 'anonymous'
    timeanddate <- as.character(Sys.time())
    analysis_name <-
      paste0(username, '_', timeanddate, '_tviblindi', packageVersion('tviblindi'), '_analysis')
  }
  
  tv <- new_tviblindi(data, labels, fcs_path, fcs_subset_idcs, analysis_name)
  
  .msg('Initialised trajectory inference analysis object ', endline = FALSE); .msg_name(analysis_name)
  .msg('Number of events: ', endline = FALSE); .msg_alt(nrow(data))
  .msg('Expression matrix dimensionality: ', endline = FALSE); .msg_alt(ncol(data)); .msg('')
  
  tv
}

new_tviblindi <- function(
  data,
  labels = NULL,
  fcs_path = NULL,
  fcs_subset_idcs = NULL,
  analysis_name = NULL
) {
  if (!is.matrix(data) || !is.numeric(data))
    stop('Data must be a numeric matrix of coordinates')
  if (length(colnames(data)) == 0)
    stop('Data must have column names set')
  if (is.null(labels)) {
    labels <- rep('ungated', nrow(data))
    warning('No labels were provided for data')
  }
  if (any(duplicated(colnames(data))))
    stop('Data column names contain duplicates')
  if (length(labels) != nrow(data))
    stop('Labels vector must have as many entries as the there are rows in the data matrix')
  if (!is.factor(labels) && !is.character(labels))
    stop('Labels must be a character or factor vector')
  
  tv                      <- new.env(hash = TRUE)
  tv$analysis_name        <- analysis_name
  tv$timestamp            <- as.character(Sys.time())
  tv$origin               <- NULL
  tv$data                 <- data
  tv$downsampling         <- NULL
  tv$denoised             <- NULL
  tv$denoising_iterations <- 0
  tv$labels               <- if (!is.list(labels)) list(default = as.factor(labels)) else labels
  tv$pseudotime           <- NULL
  tv$filtration           <- NULL
  tv$filtration_method    <- NULL
  tv$boundary             <- NULL
  tv$reduced_boundary     <- NULL
  tv$walks                <- NULL
  tv$fates                <- NULL
  tv$kNN                  <- NULL
  tv$dist                 <- NULL
  tv$trans                <- NULL
  tv$trans_kernel         <- NULL
  tv$som_xgrid            <- NULL
  tv$som_ygrid            <- NULL
  tv$clusters             <- NULL
  tv$metaclusters         <- NULL
  tv$codes                <- NULL
  tv$layout               <- NULL
  tv$layout_method        <- NULL
  tv$vae                  <- NULL
  tv$fcs_subset_idcs      <- fcs_subset_idcs
  tv$fcs                  <- fcs_path
  tv$n_walks              <- NULL
  tv$walks_kernel         <- NULL
  tv$walks_equinumerous   <- NULL
  gc()
  
  structure(tv, class = 'tviblindi')
}

.format_with_linebreaks <- function(x, items_per_line = 6) {
  m <- length(x)
  n <- m %/% items_per_line
  
  res <- ''
  
  for (idx.line in seq_len(max(n, 1))) {
    res <- paste0(res, '\n\t\t', paste(x[1:min(m, items_per_line)], collapse = ', '))
    x <- x[-(seq_len(min(m, items_per_line)))]
  }
  if (m > 0)
    res <- paste0(res, '\n\t\t', paste(x, collapse = ', '))
  res
}

print.tviblindi <- function(x) {
  .msg('tviblindi trajectory inference object')
  .msg('\tname of analysis:       ', endline = FALSE); .msg_name(x$analysis_name)
  .msg('\ttimestamp:              ', endline = FALSE); .msg_alt(x$timestamp)
  .msg('\tevents count:           ', endline = FALSE); .msg_alt(nrow(x$data))
  .msg('\tparameters count:       ', endline = FALSE); .msg_alt(ncol(x$data))
  .msg('\tdownsampling:           ', endline = FALSE); .msg_alt(if (!is.null(x$downsampling)) { x$downsampling } else { 'none' })
  
  .msg('\tlabel names:');
  for (lname in names(x$labels)) {
    .msg('\t  ', endline = FALSE); .msg_name(lname, endline = FALSE); .msg(':', endline = FALSE)
    .msg_alt(.format_with_linebreaks(levels(x$labels[[lname]])))
  }
  
  .msg('\tmarker names:           ', endline = FALSE); .msg_alt(.format_with_linebreaks(colnames(x$data)))
  .msg('\tinput FCS file:         ', endline = FALSE); .msg_alt(if (!is.null(x$fcs)) { x$fcs } else { 'none' })
  
  .msg('\tcell-of-origin indices: ', endline = FALSE)
  if (is.null(x$origin))
    .msg_alt('not set')
  else {
    .msg('')
    for (oname in names(x$origin)) {
      .msg('\t  ', endline = FALSE); .msg_name(lname, endline = FALSE); .msg(':', endline = FALSE)
      .msg_alt(x$origin[[lname]])
    }
  }
  
  .msg('\tk-NN matrix:            ', endline = FALSE); .msg_alt(if (!is.null(x$kNN)) { paste0('available, k=', ncol(x$kNN$IND) - 1) } else { 'not available' })
  .msg('\tdenoising iterations:   ', endline = FALSE); .msg_alt(x$denoising_iterations)
  .msg('\t2-d layout:             ', endline = FALSE); .msg_alt(if (!is.null(x$layout)) { paste0(x$layout_method, collapse = ', ') } else { 'not available' })
  .msg('\tdistance matrix:        ', endline = FALSE); .msg_alt(if (is.null(x$dist)) { 'not available' } else { 'available' })
  
  .msg('\ttransition models:      ', endline = FALSE)
  if (is.null(x$trans))
    .msg_alt('not constructed')
  else {
    .msg('')
    for (tname in names(x$trans)) {
      .msg('\t  ', endline = FALSE); .msg_name(tname, endline = FALSE); .msg(': ', endline = FALSE)
      .msg_alt(paste0('available, kernel=', x$trans_kernel[[tname]]))
    }
  }
  
  .msg('\tsimulated walks:        ', endline = FALSE)
  if (is.null(x$walks))
    .msg_alt('not available')
  else {
    .msg('')
    for (wname in names(x$walks)) {
      .msg('\t  ', endline = FALSE); .msg_name(wname, endline = FALSE); .msg(': ', endline = FALSE)
      .msg_alt(paste0(x$n_walks[[wname]], ' walks simulated, kernel=', x$walks_kernel[[wname]], ', equinumerous=', as.character(x$walks_equinumerous[[wname]])))
    }
  }

  .msg('\tclustering:             ', endline = FALSE); .msg_alt(if (!is.null(x$codes)) { if (x$cluster_method == 'som') paste0(x$som_xgrid, ' x ', x$som_ygrid, ' SOM grid') else paste0(x$k, ' k-means clusters') } else { 'not computed' })
  .msg('\tfiltration:             ', endline = FALSE); .msg_alt(if (!is.null(x$filtration_method)) { x$filtration_method } else { 'not available' })
}

Plot <- function(
  tv, layout_idx, labels_idx, transition_model_idx, markers, labels, pch, cex, col, plots_per_row, main, cex.main, hide_legend, ...
) UseMethod('Plot', tv)

#' Plot 2-dimensional layout of data
#'
#' Plots previously generated 2-dimensional layout of expression data in \code{tviblindi} object, with custom colour-coding of data points.
#' 
#' To colour events by pseudotime, make sure the \code{pseudotime} slot is filled first. A yellow-to-brown-to-red (increasing) colour range is used to denote pseudotime values.
#' 
#' Colouring events by labels gives a plot with a colour-coded legend.
#' Alternatively (if the previous approach results in an overly cluttered plot), he localisation of point belonging to a single label of interest can be displayed.
#'
#' Finally, yellow-to-brown-to-red heat map of expression of a single marker of interest can be displayed
#'
#' @param tv \code{tviblindi}-class object with \code{layout}
#' @param layout_idx integer or string: index of layout to use. Default value is \code{1}
#' @param labels_idx integer or string: index of per-event label vector to use. Default value is \code{1}
#' @param transition_model_idx integer or string; index of transition model to use. Default value is \code{1}
#' @param markers string or string vector: name of marker of interest for constructing a heatmap. If multiple marker names are given (maximum 32), a grid of plots is produced.  If a substring of a marker name is given and it can be matched, the most probable match will be used. Default value is \code{NULL}
#' @param labels string or string vector: name of population of interest for constructing a 2-d map of that population. If multiple population names are given (maximum 32), a grid of plots is produced. If parameter \code{markers} is not \code{NULL}, it overrides this parameter. Default value is \code{NULL}
#' @param plots_per_row integer: if a grid of multiple plots is to be generated, how many plots should fit onto a single row? Default value is \code{4}
#' @param pch integer or character: \code{pch} parameter of the \code{plot} method (plotting character or symbol). Default value is '\code{.}' (a dot)
#' @param cex numeric: \code{cex} parameter of the \code{plot} method: plotting character or symbol size. Default value is \code{3}
#' @param col string: can be '\code{labels}' or '\code{pseudotime}' to specify type of plot. '\code{labels}' is for plotting a map of all labelled cell populations; '\code{markers}' is for showing a grid of heatmaps for the expression of each marker. If either parameter \code{markers} or \code{labels} is not \code{NULL}, it overrides this parameter. Default value is '\code{labels}'.
#' @param main string: if specified, this overrides the default \code{main} parameter of \code{plot} (title of the plot, or absence thereof). Default value is \code{NULL}
#' @param cex.main number: if specified, this overrides the default \code{cex.main} parameter of \code{plot} (font size of the plot title). Default value is \code{NULL}
#'
#' @export
Plot.tviblindi <- function(
  tv,
  layout_idx = 1,
  labels_idx = 1,
  transition_model_idx = 1,
  markers = NULL,
  labels = NULL,
  pch = '.',
  cex = 3,
  col = 'labels',
  plots_per_row = 4,
  main = NULL,
  cex.main = NULL,
  hide_legend = FALSE,
  ...
) {
  if (is.null(tv$layout))
    stop('2-dimensional layout of data needs to be created first')
  if (is.numeric(layout_idx) && layout_idx > length(tv$layout))
    stop('Invalid numeric "layout_idx" value')
  if (is.character(layout_idx) && !layout_idx %in% names(tv$layout))
    stop('Invalid character "layout_idx" value')
  if (is.numeric(labels_idx) && labels_idx > length(tv$labels))
    stop('Invalid numeric "labels_idx" value')
  if (is.character(labels_idx) && !labels_idx %in% names(tv$labels))
    stop('Invalid character "labels_idx" value')
  if (!is.null(tv$pseudotime) && (is.numeric(transition_model_idx) && transition_model_idx > length(tv$pseudotime)))
    stop('Invalid numeric "transition_model_idx" value')
  if (!is.null(tv$pseudotime) && (is.character(transition_model_idx) && !labels_idx %in% names(tv$pseudotime)))
    stop('Invalid character "transition_model_idx" value')
  
  layout <- tv$layout[[layout_idx]]
  labels_per_event <- tv$labels[[labels_idx]]
  pseudotime <- tv$pseudotime[[transition_model_idx]]$res
  
  if (is.null(markers)) {
    if (!is.null(labels)) {
      for (l in labels)
        if (!l %in% levels(labels_per_event)) stop(paste0("Label '", l,"' was not found in the data.\nPlease pick one of the following labels: ", paste(as.character(levels(labels_per_event)), collapse = ', ', sep = ''), '.'))
      labels <- unique(labels)
      gating_palette <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
      gating_palette <- unlist(mapply(RColorBrewer::brewer.pal, gating_palette$maxcolors, rownames(gating_palette)))[-1]
      gating_palette <- gating_palette[1:length(labels)]
      
      xpd.old <- par()$xpd
      mar.old <- par()$mar
      par(xpd = TRUE, mar = c(2, 2, 2, 2))
      
      .msg('Plotting select labelled populations')
      
      plot(
        layout,
        col = scales::alpha('black', .3),
        axes = FALSE,
        xlab = '',
        ylab = '',
        pch = pch,
        cex = cex,
        main = if (is.null(main)) { '' } else { main },
        cex.main = if (is.null(cex.main)) { 2 } else { cex.main }
      )
      
      for (i in seq_along(labels)) {
        label <- labels[i]
        points(
          layout[labels_per_event == label, 1],
          layout[labels_per_event == label, 2],
          col = scales::alpha(gating_palette[i], .5),
          pch = pch,
          cex = cex + 0.5)
      }

      legend('bottomleft', legend = labels, fill = gating_palette, bty = 'n', cex = .8, xpd = TRUE)      
      
      par(xpd = xpd.old, mar = mar.old)
      
    } else {
      if (col[1] == 'pseudotime') {
        if (is.null(tv$pseudotime))
          stop('Pseudotime needs to be computed first')
        psc <- as.numeric(as.factor(pseudotime))
        psc <- psc / max(psc)
        psc <- psc * 10000 + 1
        col <- gplots::colorpanel(10500, low = 'yellow', mid = 'brown', high = 'red')
        
        .msg('Plotting events coloured by pseudotime values')
        
        plot(
          layout,
          col = scales::alpha(col[psc], 0.2),
          pch = pch,
          cex = cex,
          axes = FALSE,
          xlab = '',
          ylab = '',
          main = if (is.null(main)) { '' } else { main },
          cex.main = if (is.null(cex.main)) { 1.5 } else { cex.main }
        )
        
      } else if (col[1] == 'labels') {
        
        ## Set up colour palette for plotting annotated populations
        gating_palette <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
        gating_palette <- unlist(mapply(RColorBrewer::brewer.pal, gating_palette$maxcolors, rownames(gating_palette)))[-c(1, 12, 16)]
        gating_palette <- gating_palette[1:length(unique(labels_per_event))]
        gating_colour_vector <- gating_palette[as.numeric(as.factor(labels_per_event))]
        colours.aligned <- gating_palette[1:length(unique(labels_per_event))]
        
        xpd.old <- par()$xpd
        mar.old <- par()$mar
        par(xpd = TRUE, mar = c(2, 2, 2, 2))
        
        message('Plotting all labelled populations')
        
        # Plot ungated events distinctly
        labels.aligned <- sort(unique(labels_per_event))
        which.ungated <- which(labels.aligned %in% c('ungated', '*ungated*'))
        if (length(which.ungated) == 1) {
          label.ungated <- labels.aligned[which.ungated]
          idcs.ungated <- which(labels_per_event == label.ungated)
          colours.aligned[which.ungated] <- gating_colour_vector[idcs.ungated] <- 'black'
          plot(layout[idcs.ungated, ], col = scales::alpha('black', .4), axes = FALSE, xlab = '', ylab = '', pch = pch, cex = cex, main = if (is.null(main)) { '' } else { main }, cex.main = if (is.null(cex.main)) { 1.5 } else { cex.main })
          points(layout[-idcs.ungated, 1], layout[-idcs.ungated, 2], col = scales::alpha(gating_colour_vector[-idcs.ungated], 0.4), pch = pch, cex = cex)
        } else {
          plot(
            layout,
            col = scales::alpha(gating_colour_vector, .35),
            axes = FALSE,
            xlab = '',
            ylab = '',
            pch = pch,
            cex = cex,
            main = if (is.null(main)) { '' } else { main },
            cex.main = if (is.null(cex.main)) { 1.5 } else { cex.main }
          )
        }
        
        if (!hide_legend)
          legend('bottomleft', legend = labels.aligned, fill = colours.aligned, bty = 'n', cex = .8, xpd = TRUE)
        
        par(xpd = xpd.old, mar = mar.old)
      }
    }
  } else {
    par(xpd = TRUE, mar = c(6, 3, 3, 3))
    for (i in 1:length(markers)) {
      m <- markers[i]
      
      .msg(paste0('Plotting expressions of marker ', m))
      
      M <- colnames(tv$data)
      if (!m %in% M) {
        j <- grep(m, M, fixed = TRUE)
        if (length(j) == 0)
          stop(paste0("Marker '", m,"' was not found in the data.\nPlease pick one of the following markers: ", paste(M, collapse = ', ', sep = ''), '.'))
        else {
          if (length(j) > 1) {
            # Preferentially pick the match which either begins or ends with the substring (eg. if 'CD4' is the substring, match it to '145Nd_CD4' not '89Y_CD45').
            which.left <- which(substr(M, 1, nchar(m)) == m)
            which.right <- which(substr(M, nchar(M) - nchar(m) + 1, nchar(M)) == m)
            if (length(which.left  == 1)) j <-  which.left
            if (length(which.right == 1)) j <- which.right
            if (length(which.left   > 1)) j <-  which.left[1]
            if (length(which.right  > 1)) j <- which.right[1]
            if (length(j)           > 1)  j <-           j[1]
          }
          message(paste0("Marker '", m, "' was not found in the data. It was matched to '", M[j], "'"))
          markers[i] <- M[j]
        }
      }
    }
    maxval <- max(tv$data)
    minval <- min(tv$data)
    col <- gplots::colorpanel(10500, low = 'yellow', mid = 'brown', high = 'red')
    
    mfrow.old <- par()$mfrow
    n.markers <- length(markers)
    ppr <- plots_per_row
    if (n.markers <= ppr)
      par(mfrow = c(1, n.markers))
    else
      par(mfrow = c((n.markers - 1) %/% ppr + 1, ppr))
    
    for (marker in markers) {
      n <- 10500
      inc <- (maxval - minval) / (n - 1)
      ticks <- minval + (0:(n - 1)) * inc
      
      col_idcs <- cut(tv$data[, marker], breaks = ticks, labels = FALSE)
      
      stats <- summary(tv$data[, marker])[-4]
      sub <- paste(c('min', '1st quartile', 'median', '3rd quartile', 'max'), round(stats, 3), sep = ': ', collapse = '\n')
      plot(layout, col = scales::alpha(col[col_idcs], 0.2), pch = pch, cex = cex, axes = FALSE, xlab = '', ylab = '', main = if (is.null(main)) { marker } else { main }, cex.main = if (is.null(cex.main)) { 2 } else { cex.main })
      box(which = 'plot', lty = 'solid', col = '#e6e6e6', lwd = 1.5)
      adj.old <- par()$adj
      par(adj = 1)
      title(sub = sub, cex.sub = 1.6)
      par(adj = adj.old)
    }
    
    par(mfrow = mfrow.old)
  }
}

#' Add a vector of labels per event to a \code{tviblindi} object
#'
#' Adds a vector of labels per each row of the \code{data} matrix of a \code{tviblindi} object.
#' (Multiple alternative label vectors can be supplied, and multiple TI analyses can be conducted, with separate sets of simulated random walks for each labelling of the data.)
#'
#' The default vector of labels is named '\code{default}'. Use '\code{default}' as a name for your new labels vector to overwrite this.
#'
#' @param tv \code{tviblindi}-class object
#' @param labels string or factor vector: label for each event (row in \code{tv$data})
#' @param name string: name of the labels vector
#'
#' @export
AddLabelsVector.tviblindi <- function(
  tv,
  labels,
  name = 'default'
) {
  if (length(labels) != nrow(tv$data))
    stop('Invalid length of label vector')
  if (!is.null(tv$labels[[name]])) {
    ans <- readline(prompt = paste0('Label vector "', name, '" already exists. Overwrite? (Y/n) '))
    if (ans != 'Y')
      return(NULL)
  }
  tv$labels[[name]] <- as.factor(labels)
  invisible(tv)
}

AddLabelsVector <- function(tv, labels, name) UseMethod('AddLabelsVector', tv)

#' Make a deep copy of a \code{tviblindi}
#'
#' Adds a vector of labels per each row of the \code{data} matrix of a \code{tviblindi} object.
#' (Multiple alternative label vectors can be supplied, and multiple TI analyses can be conducted, with separate sets of simulated random walks for each labelling of the data.)
#'
#' The default vector of labels is named '\code{default}'. Use '\code{default}' as a name for your new labels vector to overwrite this.
#'
#' @param tv \code{tviblindi}-class object
#' @param labels string or factor vector: label for each event (row in \code{tv$data})
#' @param name string: name of the labels vector
#'
#' @export
Copy.tviblindi <- function(tv) {
  res <- new.env()
  for (obj in ls(tv))
    res[[obj]] <- if(is.environment(obj)) env.copy(tv[[obj]], all.names) else tv[[obj]]
  structure(res, class = 'tviblindi')
}

Copy <- function(tv) UseMethod('Copy', tv)