## tviblindi
## Shiny interface module
## Server function

require(shinyBS)

#### Shiny interface for tviblindi: server function

shiny_server <- function(input, output, session) {
  
  cancel.onSessionEnded <- session$onSessionEnded(function() {
    message('tviblindi Shiny session ended')
    stopApp()
  })
  message(paste0('Running tviblindi Shiny UI, working directory is ', getwd()))
  
  ## Load up tviblindi analysis S3 object
  INPUTS_DIRECTORY <- 'tviblindi_tmp'
  tv_name          <- readRDS(file.path(INPUTS_DIRECTORY, 'tv.RDS'))
  tv               <- get(tv_name, parent.env(environment()))
  
  input_ff         <- if (!is.null(tv$fcs)) { flowCore::read.FCS(tv$fcs) } else { make_valid_fcs(exprs = tv$data) }
  layout           <- if (class(tv$layout) == 'list') { tv$layout } else { list(default = tv$layout) }
  labels           <- if (class(tv$labels) == 'list') { tv$labels } else { list(default = tv$labels) }
  event_sel        <- tv$events_sel # selected events from input FCS file
  markers          <- colnames(tv$data)
  
  if ((is.null(event_sel) && nrow(input_ff) != nrow(tv$data)) || (!is.null(event_sel) && length(event_sel) != nrow(tv$data))) {
    stop('Number of events in expression matrix incongruent with dimensionality of input FCS file. Did you misuse the event_sel parameter?')
  }
  
  ## Reactive values
  react                              <- reactiveValues()
  react$pseudotime                   <- tv$pseudotime
  react$persistence                  <- NULL
  react$persistence_diagram          <- NULL
  react$representations              <- NULL # trajectories representations
  react$representations.reduced      <- NULL # trajectory representations only using marked homology classes
  # Dimred layout
  react$layout_name                  <- names(layout)[1]
  react$layout_pointsize             <- .065
  # Labels vector
  react$labels_name                  <- names(labels)[1]
  # Terminal nodes selection
  react$termini_selected             <- NULL # selected terminal nodes idcs
  react$termini_marked               <- NULL # marked terminal nodes idcs
  react$termini_processed            <- NULL
  # Homology classes selection
  react$persistence_selected         <- NULL # selected homology classes idcs
  react$persistence_marked           <- NULL # marked homology classes idcs
  react$persistence.death_birth_ratio <- FALSE
  # Trajectories dendrogram
  react$dendrogram                   <- NA   # dendrogram for clustering of trajectories by homology classes
  react$dendrogram_zoom              <- NA
  react$dendrogram_zoomed            <- FALSE
  react$dendrogram_data              <- NA
  react$dendrogram_zoom_data         <- NA
  react$dendrogram_ready             <- FALSE
  react$dendrogram_zoom_ready        <- FALSE
  react$dendrogram_selection         <- NULL
  react$dendrogram_selected_leaves   <- NULL
  react$dendrogram_selected_idcs     <- NULL
  react$dendrogram_marked_leaves.A   <- NULL
  react$dendrogram_zoom_idcs         <- NULL
  react$dendrogram_redraw_highlights <- FALSE
  react$dendrogram_redraw_zoom       <- FALSE
  react$dendrogram_zoom_redraw_highlights <- FALSE
  react$dendrogram_leaf_perc_cutoff      <- NULL
  react$dendrogram_classes           <- list()
  react$dendrogram_zoom_classes      <- list()
  react$dendrogram_labels            <- NA
  react$dendrogram_zoom_labels       <- NA
  react$dendrogram_zoom_active       <- FALSE
  react$dendrogram_draw_zoomed_in    <- FALSE
  # 2D trajectories layout
  react$layout_trajectories_flip_colours            <- FALSE
  react$layout_trajectories_highlight_in_background <- FALSE
  # Marker expression and population trackers
  react$trackers_n_segments          <- 20 # number of pseudotime segments in markers tracking
  react$trackers_scaling_exponent    <- 1 # scaling exponent for pseudotime segmentation in markers tracking
  react$trackers_large_base_size     <- FALSE
  react$tracked_markers.A            <- NULL # names of markers of interest in group A
  react$tracked_markers_stats.A      <- NULL
  react$tracked_markers_stats.B      <- NULL
  react$tracked_markers_parts_to_remove.A <- NULL
  react$tracked_markers_parts_to_remove.B <- NULL
  react$tracked_markers_ready.A      <- FALSE
  react$tracked_markers_ready.B      <- FALSE
  react$tracked_populations_ready.A      <- FALSE
  react$tracked_populations_ready.B      <- FALSE
  react$tracked_markers_pseudotime_bounds.A <- NULL
  react$tracked_markers_pseudotime_bounds.B <- NULL
  react$tracked_markers_last_removed.A <- NULL
  react$tracked_markers_last_removed.B <- NULL
  react$pseudotime_highlight_bounds  <- NULL
  react$selected_trajectory_points <- NULL
  react$tracked_populations.A        <- NULL
  react$tracked_populations.B        <- NULL
  react$tracked_populations_stats.A  <- NULL
  react$tracked_populations_stats.B  <- NULL
  react$tracked_populations_log2_transform <- FALSE
  # Marked trajectories manipulation & export
  react$trajectories_marked.A        <- NULL
  react$dendrogram_marked_leaves.A <- NULL
  react$trajectories_marked.B        <- NULL
  react$dendrogram_marked_leaves.B <- NULL
  react$trajectories_junk.A          <- integer(0)
  react$trajectories_junk.B          <- integer(0)
  react$trajectories_group           <- 'A'
  react$output_ff                    <- NULL
  react$trajectories_pinned_batches_count <- 0
  react$trajectories_random_walks    <- NULL
  react$trajectories_to_pin          <- NULL
  react$trajectories_pinned          <- NULL
  # Images export
  react$image_export.termini           <- FALSE
  react$image_export.persistence       <- FALSE
  react$image_export.dendrogram        <- FALSE
  react$image_export.dendrogram_zoom   <- FALSE
  react$image_export.trajectories      <- FALSE
  react$image_export.tracked_markers.A <- FALSE
  react$image_export.tracked_markers.B <- FALSE
  react$image_export.tracked_populations.A <- FALSE
  react$image_export.tracked_populations.B <- FALSE
  react$image_export_format            <- 'PNG'
  
  ## Scale 2D projection
  for (idx.layout in names(layout)) {
    layout[[idx.layout]][, 1] <- layout[[idx.layout]][, 1] - min(layout[[idx.layout]][, 1]); layout[[idx.layout]][, 1] <- layout[[idx.layout]][, 1] / max(layout[[idx.layout]][, 1])
    layout[[idx.layout]][, 2] <- layout[[idx.layout]][, 2] - min(layout[[idx.layout]][, 2]); layout[[idx.layout]][, 2] <- layout[[idx.layout]][, 2] / max(layout[[idx.layout]][, 2])
  }
  
  # Compute pseudotime colouring
  psc   <- as.numeric(as.factor(tv$pseudotime$res))
  psc   <- psc / max(psc)
  psc   <- psc * 10000 + 1
  pal.hex  <- gplots::colorpanel(10500, low = 'yellow', mid = 'brown', high = 'red')
  pal.rgba <- col2rgb(pal.hex, alpha = 0.7)
  cols <- pal.hex[psc]
  cols_as_matrix <- pal.rgba[, psc]
  
  layout.df <- lapply(layout, function(l) data.frame(l) %>% setNames(c('X', 'Y'))) # for compatibility with brushOpts
  names(layout.df) <- names(layout)
  
  gating_palette <- c(RColorBrewer::brewer.pal(8, 'Dark2'), RColorBrewer::brewer.pal(12, 'Paired')[-11], RColorBrewer::brewer.pal(9, 'Set1')[-6], RColorBrewer::brewer.pal(8, 'Accent')[5:8])
  labels.unique <- lapply(labels, unique) %>% setNames(names(labels))
  labels.aligned <-
    colours.aligned <-
    gating_colour_vectors <- vector(mode = 'list', length = length(labels)) %>% setNames(names(labels))
  for (idx.labels in 1:length(labels)) {
    # Re-order annotated population labels by pseudotime
    average_pseudotime <- c()
    for (label in labels.unique[[idx.labels]]) {
      idcs        <- which(labels[[idx.labels]] == label)
      idcs.sample <- sample(idcs, min(100000, length(idcs)))
      average_pseudotime <- c(average_pseudotime, mean(tv$pseudotime$res[idcs.sample]))
    }
    label_levels_order <- order(average_pseudotime)
    labels.aligned[[idx.labels]] <- labels.unique[[idx.labels]][label_levels_order]
    
    gating_colour_vectors[[idx.labels]] <- gating_palette[1:length(labels.unique[[idx.labels]])][as.numeric(as.factor(labels[[idx.labels]]))]
    colours.aligned[[idx.labels]] <- sapply(labels.aligned[[idx.labels]], function(l) gating_colour_vectors[[idx.labels]][which(labels[[idx.labels]] == l)[1]])
  }
  
  observeEvent(input$input_labels_name, {
    react$labels_name <- input$input_labels_name
  })
  
  ## Find all trajectories' terminal nodes
  termini              <- c(tv$walks$v[tv$walks$starts[-1] - 1], tail(tv$walks$v, 1))
  termini.unique       <- unique(termini)
  
  ### OUTPUTS & OBSERVERS
  
  ## Reset all selectors if window is resized (otherwise this causes bugs)
  observeEvent(input$windowSizeChange, {
    session$resetBrush('selector_termini')
    session$resetBrush('selector_persistence')
    session$resetBrush('selector_dendrogram')
    session$resetBrush('selector_dendrogram_zoom')
    session$resetBrush('selector_tracked_markers.A')
    session$resetBrush('selector_tracked_markers.B')
  })
  
  ## Display dimension reduction methods and labels vector names in dropdown menu
  observe({
    updateSelectInput(session, 'input_dimred_method', choices = names(layout))
    updateSelectInput(session, 'input_labels_name', choices = names(labels))
  })
  observeEvent(input$input_dimred_method, {
    react$layout_name <- input$input_dimred_method
  })
  
  observeEvent(input$btn_layout_pointsize, {
    react$layout_pointsize <- as.numeric(input$btn_layout_pointsize)
  })
  
  ## Display FCS file name
  if (!is.null(tv$analysis_name)) {
    output$text_analysis_name <- renderText(if (tv$analysis_name != FALSE) { tv$analysis_name } else { 'Untitled tviblindi analysis' })
  }
  if (!is.null(tv$fcs)) {
    output$text_fcs_name <- renderText(paste0('Linked to FCS file: ', tv$fcs))
  }
  observeEvent(input$btn_help, {
    showModal(modalDialog(
      title = HTML('<h3>How to use the <i>tviblindi</i> Shiny interface</h3>'),
      size = 'l',
      HTML("This interface is a tool to allow for discrimination within a set of canonical developmental trajectories in input data, as well as between canonical and aberrant trajectories. Some trajectories will be biologically relevant, others not.
             <br><br>
             First, select terminal nodes of simulated random walks in the left pane (tab <i>Terminal nodes selection</i>). This is done by clicking and drawing a rectangular selection. To mark selected points, press the <b>PLUS</b> button underneath. To remove marked points (and start selecting over again), press the <b>FIRE</b> button underneath. Once you are satisfied with your selection (the marked terminal nodes will be clustered together!), press the <b>THUMBS-UP</b> button to compute the relevant triangulation and create trajectory representations. Additionally, press the <b>ENLARGE</b> button in the right part of the left panel to view the dimension reduction layout coloured according to annotated cell populations.
             <br><br>
             Second, select the tab <i>Homology classes by persistence selection</i> in the left pane. Here, select relevant homology classes to be used for hierarchical classification of trajectories. (For our analysis, these correspond roughly to more or less significant holes in the high-dimensional data.) Persistence increases upward. Again, you can make multiple piecewise selections by consecutively drawing rectangles and pressing <b>PLUS</b> to mark the selected points, or you can discard the marked points by pressing the <b>FIRE</b> button.
             <br><br>
             Third, select and mark trajectories of interest using the dendrogram which appears in the middle pane. Each bifurcation in the dendrogram corresponds to difference in path classification with regard to a single homology class. The horizontal coordinate of each bifurcation corresponds to filtration value of addition of a simplex associated with the death of that homology class during filtration. In other words, as we move toward the right, trajectories as clustered together based on the differences in how they circumnavigate around increasingly prominent holes in the high-dimensional space. During the marking of trajectories using the <b>PLUS</b> button, be sure to select the desired group (<b>A</b> or <b>B</b>) for distinguishing between two collections of marked trajectories.
             Upon adding any branches of the dendrogram to group A, a blue rectangle is drawn over the dendrogram, showing the marked branches. A red rectangle is displayed for branches in group B.
             If the dendrogram is particularly complex to navigate, select an area of interest and click the <b>MAGNIFYING GLASS</b> button to generate a zoomed-in view of that area. A green rectangle will appear on the right border of the plot, showing which area was magnified.
             Then, switch from the tab <i>Whole dendrogram</i> to the tab <i>Zoom</i> to interact with the zoomed-in area. To change the magnified area, return to the <i>Whole dendrogram</i> tab, change your selection and click the <b>MAGNIFYING GLASS</b> button again. To clear the magnified selection, leave your selection empty and then click the <b>MAGNIFYING GLASS</b>.
             <br><br>
             Fourth, inspect the projection of marked trajectories in the 2D layout in the right pane. Category <b>A</b> trajectories are drawn in blue, whereas category <b>B</b> trajectories are drawn in blue. By default, red trajectories are drawn on the top. To flip this ordering, press the <b>FLIP</b> button (with the black-white circular icon) beneath the 2D layout plot.
             <br><br>
             Fifth, inspect the progression of marker expression in either category of trajectories by selecting the <i>Marker expression tracking</i> tab and entering marker(s) of interest in the text boxes below the 2D projection (category <b>A</b> is above category <b>B</b> here). By default, we separate the progression into 20 equally large segments by exponentially scaled pseudotime values. Both the scaling exponent and the number of segments can be adjusted (however, the default settings should give sensible results in most cases). If a single marker of interest is entered, all the trajectories are included in the diagram. In addition, you can remove (un-mark) some trajectories by drawing a selection rectangle in the progression diagram and pressing the <b>LIGHTNING BOLT</b> button. All trajectories passing through the drawn area are 'zapped' (the number of 'zapped' trajectories is printed in the bottom part of the center pane). You cannot zap all pathways in a category. You can undo a single last zap by pressing the <b>BACKWARD</b> button.
             If multiple markers of interest are chosen, the mean value for each segment for each marker is displayed. In that case zapping is not possible. If you are interested in certain pseudotime segments (for example, corresponding to some sudden spike in marker expressions), select points in those segments and press the <b>FLAG</b> button to temporaritly overlay the trajectories layout with a highlight of those points which belong in the corresponding psedotime segments. After pressing the <b>FLAG</b> button, the selection is automatically cleared, therefore clicking the button again will undo the highlighting. Alternatively, press the <b>CROSS</b> button under the 2D trajectories layout to remove the highlighting. In order to switch between a background and foreground positioning of the highlight, press the <b>GLASSES</b> button. Pseudotime segments highlighting works both when viewing single and multiple markers progressions.
             <br><br>
             Sixth, inspect the counts of select annotated cell populations in each pseudotime segment by selecting the <i>Population tracking</i> tab. This is analogical to the marker tracking part. If you wish to display log2 of cell counts, check the <i>log2</i> checkbox.
             <br><br>
             Seventh, export an enhanced FCS file. In the simplest use-case, you can press the <b>SAVE</b> button in the middle pane straightaway (without following any of the above steps) to append two artificial channels to your FCS, containing the 2D layout coordinates displayed in the left and right panes. If you want to append marked trajectories in either of the categories, along with pseudotime values, press the <b>PIN</b> button next to the header for either category. Then, press the <b>SAVE</b> button. To un-pin all the batches of vertices pinned so far, press the <b>GARBAGE</b> button next to it.
             <br><br>
             Eighth, use the <b>PICTURE</b> buttons under the pseuodotime layout, modified persistence diagram, dendrogram, trajectory layout or tracked marker expression diagrams to generate PNG or SVG images in the working directory. This way, you can document your analysis and create a report.
             To choose whether to export a PNG or an SVG file, use the switch at the bottom of the left pane. PNG files are raster images. These will typically save space on your drive, but are resolution-dependent (they appear pixelated when enlarged). SVG files, on the other hand, typically take up more space (for large single-cell data layouts, possibly hunreds of megabytes), but are resolution-independent. This is especially suitable for further image processing."
      )))
  })
  
  ## Terminal nodes picker
  # Layout
  output$plot_termini <- renderPlot({
    
    if (react$image_export.termini) {
      if (react$image_export_format == 'SVG') {
        svg(filename = paste0('Termini_', Sys.time(), '.svg'))
      } else {
        png(filename = paste0('Termini_', Sys.time(), '.png'), width = 1000, height = 900)
      }
      par(mar = c(0, 0, 0, 0))
      layout.raster <- scattermore(
          xy = layout[[react$layout_name]],
          cex = react$layout_pointsize,
          rgba = cols_as_matrix,
          xlim = c(0, 1.2),
          ylim = c(0, 1)
      )
      plot(layout.raster)
      layout.raster.origin <- scattermore(
        xy = layout[[react$layout_name]][tv$origin, , drop = FALSE],
        cex = if (react$image_export.termini && react$image_export_format != 'SVG') { 5 } else { 3 },
        rgba = col2rgb('purple', alpha = .75),
        xlim = c(0, 1.2),
        ylim = c(0, 1)
      )
      layout.raster.termini <- scattermore(
        xy = layout[[react$layout_name]][termini.unique, , drop = FALSE],
        cex = if (react$image_export.termini && react$image_export_format != 'SVG') { 5 } else { 3 },
        rgba = col2rgb('grey', alpha = .75),
        xlim = c(0, 1.2),
        ylim = c(0, 1)
      )
      plot(layout.raster.origin, add = TRUE)
      plot(layout.raster.termini, add = TRUE)
      if (tv$ShowAllFates) {
        layout.raster.fates <- scattermore(
          xy = layout[[react$layout_name]][unique(tv$states), , drop = FALSE],
          cex = if (react$image_export.termini && react$image_export_format != 'SVG') { 5 } else { 3 },
          rgba = col2rgb('grey', alpha = .75),
          xlim = c(0, 1.2),
          ylim = c(0, 1)
        )
        plot(layout.raster.fates, add = TRUE)
      }
    } else {
      par(mar = c(1, 1, 1, 1))
      plot(
        layout[[react$layout_name]],
        col = scales::alpha(cols, .7),
        axes = FALSE, xlab = '', ylab = '',
        pch = if (react$layout_pointsize < 0.5) { 16 } else { 20 },
        cex = if (react$image_export.termini && react$image_export_format != 'SVG') { react$layout_pointsize + .3 } else { react$layout_pointsize },
        xlim = c(0, 1.2), ylim = c(0, 1)
      )
      # Plot origin and terminal nodes
      points(
        layout[[react$layout_name]][tv$origin, , drop = FALSE],
        col = scales::alpha('purple', .75),
        cex = if (react$image_export.termini && react$image_export_format != 'SVG') { 5 } else { 3 },
        pch = 15
      )
      if (!is.null(tv$ShowAllFates) && tv$ShowAllFates) {
        points(
          layout[[react$layout_name]][unique(tv$states), , drop = FALSE],
          col = scales::alpha('grey', .75),
          cex = if (react$image_export.termini && react$image_export_format != 'SVG') { 3 } else { 1.5 },
          pch = 20
        )
      }
      points(
        layout[[react$layout_name]][termini.unique, , drop = FALSE],
        col = scales::alpha('black', .5),
        cex = if (react$image_export.termini && react$image_export_format != 'SVG') { 5 } else { 3 },
        pch = 20
      )
    }
    
    rasterImage(as.raster(matrix(colorRampPalette(c('yellow', 'brown', 'red'))(20), nc = 1)), 1.07, 0, 1.1, 0.9)
    text(x = 1.15, y = seq(0.01, 0.89, l = 3), labels = c('late', 'mid', 'early'))
    text(x = 1.10, y = 0.98, labels = 'PSEUDOTIME', font = 2)
    
    if (react$image_export.termini) {
      dev.off()
      react$image_export.termini <- FALSE
    }
  })
  
  # Buttons
  observeEvent(input$btn_termini_clear_termini, {
    react$termini_marked          <- data.frame()
    output$log_persistence_marked <- renderText('')
  })
  
  observeEvent(input$btn_termini_update_walks_by_termini, {
    react$termini_processed <- react$termini_marked
    if (length(react$termini_marked) > 0) {
      updated <- .update_walks_by_termini(tv                = tv,
                                          pseudotime        = react$pseudotime,
                                          marked_termini    = react$termini_processed,
                                          termini_per_path  = termini,
                                          death_birth_ratio = react$persistence.death_birth_ratio)
      react$trajectories_random_walks <- updated$random_walks
      react$representations           <- updated$repre
      react$persistence               <- updated$pers
      react$persistence_diagram       <- updated$pers_diag
      react$dendrogram_ready          <- TRUE
      react$persistence_selected      <- NULL
      react$persistence_marked        <- NULL
      output$log_persistence_marked   <- renderPrint({ cat('Freshly marked terminal nodes. No homology classes marked yet\n') })
      react$trajectories_marked.A     <- NULL
      react$trajectories_marked.B     <- NULL
      react$tracked_markers_last_removed.A <- NULL
      react$tracked_markers_last_removed.B <- NULL
    }
  })
  
  observeEvent(input$btn_termini_export_image, {
    react$image_export.termini <- TRUE
  })
  
  # Logs
  output$log_termini_selected <- renderPrint({
    react$termini_selected <- as.numeric(rownames(brushedPoints(layout.df[[react$layout_name]][termini.unique, ], input$selector_termini, xvar = 'X', yvar = 'Y')))
    cat(sapply(react$termini_selected, function(s) paste0(s, ' (label: ', labels[[react$labels_name]][s], '; ', sum(termini == s), ' walks terminate here)')), sep = '\n')
  })
  observeEvent(input$btn_termini_mark_termini, {
    react$termini_marked <- unique(unlist(c(react$termini_marked, react$termini_selected)))
    output$log_termini_marked <- renderPrint({ cat(sapply(react$termini_marked, function(s) paste0(s, ' (label: ', labels[[react$labels_name]][s], '; ', sum(termini == s), ' walks terminate here)')), sep = '\n') })
  })
  
  ## Gating layout pop-up
  output$plot_gating_layout <- output$plot_gating_termini_layout <- renderPlot({
    gating_layout.colours <- gating_colour_vectors[[react$labels_name]]
    pch <- rep(20, length(labels[[react$labels_name]])) # symbols (can be set to different for each population here)
    par(xpd = TRUE, mar = c(2, 2, 2, 50))
    layout.raster <- scattermore(
      xy = layout[[react$layout_name]],
      cex = rep(react$layout_pointsize * 2 + .8, nrow(layout[[react$layout_name]])),
      rgba = col2rgb(gating_colour_vectors[[react$labels_name]], alpha = 0.4),
      size = c(900, 900),
      xlim = c(0, 1),
      ylim = c(0, 1)
    )
    plot(layout.raster)
    
    layout.raster.origin <- scattermore(
      layout[[react$layout_name]][tv$origin, , drop = FALSE],
      rgba = col2rgb('#ded00b', alpha = 1),
      cex = 8,
      size = c(900, 900),
      xlim = c(0, 1),
      ylim = c(0, 1)
    )
    layout.raster.termini <- scattermore(
      layout[[react$layout_name]][termini.unique, , drop = FALSE],
      rgba = col2rgb(rep('#525252', length(termini.unique)), alpha = 1),
      cex = 8,
      size = c(900, 900),
      xlim = c(0, 1),
      ylim = c(0, 1)
    )
    suppressWarnings(plot(layout.raster.origin, add = TRUE))
    suppressWarnings(plot(layout.raster.termini, add = TRUE))
    par(mfrow = c(1, 1))
    suppressWarnings(legend(
      x = 965, y = 825,
      legend = labels.aligned[[react$labels_name]],
      fill = colours.aligned[[react$labels_name]],
      cex = 1.5,
      bty = 'n'
    ))
  })
  
  ## Homology class picker (persistence)
  # Layout
  output$plot_persistence <- renderPlot({
    if (!is.null(react$persistence)) {
      quantiles <- quantile(react$persistence_diagram$yplot)
      sizes <- max(1, sapply(react$persistence_diagram$yplot, function(val) first(which(quantiles >= val)) - 1))
      if (react$image_export.persistence) {
        if (react$image_export_format == 'SVG') {
          svg(filename = paste0('Persistence_', Sys.time(), '.svg'))
        } else {
          png(filename = paste0('Persistence_', Sys.time(), '.png'), width = 700, height = 700)
        }
        g <- ggplot(react$persistence_diagram, aes(x = xplot, y = yplot, colour = -yplot, size = yplot)) +
          scale_size_continuous(range = c(.2, 8)) +
          scale_colour_gradientn(colours = rainbow(5)) +
          geom_point(size = if (react$image_export_format == 'SVG') { 3.5 + sizes_extra * 1 } else { 5.0  + sizes_extra * 1 }, alpha = .7) +
          theme_light() + theme(legend.position = 'none') +
          xlab('(Birth + Death) / 2')
        if (react$persistence.death_birth_ratio) {
          g <- g + ylab('Death / Birth') 
        } else {
          g <- g + ylab('(Death - Birth) / 2')
        }
        if (react$image_export_format != 'SVG') {
          g <- g + theme(axis.text=element_text(size = 16),
                         axis.title=element_text(size = 20))
        }
        plot(g)
        dev.off()
        react$image_export.persistence <- FALSE
      }
      
      g <- ggplot(react$persistence_diagram, aes(x = xplot, y = yplot, colour = -yplot, size = yplot)) +
        scale_colour_gradientn(colours = rainbow(5)) +
        scale_size_continuous(range = c(.2, 8)) +
        geom_point(stroke = 1, alpha = .7) +
        theme(legend.position = 'none',
              panel.background = element_rect(fill = '#f2f2f2',
                                              colour = '#f2f2f2',
                                              size = 0.5, linetype = 'solid')) +
        xlab('(Birth + Death) / 2')
      if (react$persistence.death_birth_ratio) {
        g <- g + ylab('Death / Birth') 
      } else {
        g <- g + ylab('(Death - Birth) / 2')
      }
      g
    }
  })
  # Buttons & logs
  output$log_persistence_selected <- renderPrint({
    if (!is.null(react$persistence)) {
      react$persistence_selected <- brushedPoints(react$persistence_diagram, input$selector_persistence, xvar = 'xplot', yvar = 'yplot')
      graphical_columns          <- which(colnames(react$persistence_selected) %in% c('xplot', 'yplot'))
      if (nrow(react$persistence_selected) < 16) {
        print(react$persistence_selected[, -graphical_columns])
      } else {
        print(react$persistence_selected[1:15, -graphical_columns])
        cat('...\n')
      }
    }
  })
  
  observeEvent(input$btn_persistence_mark_classes, {
    graphical_columns        <- which(colnames(react$persistence_selected) %in% c('xplot', 'yplot'))
    react$persistence_marked <- .row_unique(rbind(react$persistence_marked, react$persistence_selected[, -graphical_columns]))
    react$persistence_marked <- react$persistence_marked[order(as.numeric(rownames(react$persistence_marked))), ]
    output$log_persistence_marked <- renderPrint({
      if (nrow(react$persistence_marked) < 16) {
        print(react$persistence_marked)
      } else {
        print(react$persistence_marked[1:15, ])
        cat('...\n')
      }
    })
    react$trajectories_marked.A <- NULL
    react$trajectories_marked.B <- NULL
    session$resetBrush('selector_persistence')
    react$dendrogram_ready <- (!is.null(react$trajectories_random_walks) && !is.null(react$persistence_marked))
  })
  observeEvent(input$btn_persistence_clear_classes, {
    react$dendrogram            <- NA
    react$dendrogram_zoom       <- NA
    react$dendrogram_ready      <- FALSE
    react$dendrogram_zoom_ready <- FALSE
    react$persistence_marked      <- data.frame()
    output$log_persistence_marked <- renderText('No homology classes marked')
    react$trajectories_marked.A   <- NULL
    react$trajectories_marked.B   <- NULL
    react$dendrogram_zoom_idcs    <- NULL
    session$resetBrush('selector_dendrogram')
    session$resetBrush('selector_dendrogram_zoom')
  })
  observeEvent(input$btn_persistence_export_image, {
    react$image_export.persistence <- TRUE
  })
  
  observeEvent(input$switch_persistence_ratio, {
    react$persistence.death_birth_ratio <- input$switch_persistence_ratio
    react$persistence_diagram       <- .persistence_diagram(react$persistence, react$persistence.death_birth_ratio)
  })
  
  ## Image export format switch
  
  observeEvent(input$btn_image_export_format, {
    react$image_export_format <- input$btn_image_export_format
  })
  
  
  ## Trajectories dendrogram
  # Plot
  output$plot_dendrogram <- renderPlot({
    session$resetBrush('selector_dendrogram_zoom')
    react$dendrogram_zoom_active <- FALSE
    par(mar = c(0, 0, 0, 0))
    perc <- react$dendrogram_leaf_perc_cutoff
    dendrogram_plotted   <- FALSE
    dendrogram_plot      <- NULL
    highlights.A         <- react$dendrogram_marked_leaves.A
    highlights.B         <- react$dendrogram_marked_leaves.B
    dendrogram_zoom_idcs <- react$dendrogram_zoom_idcs
    if (react$dendrogram_ready) {
      idcs_selected <- as.numeric(unlist(rownames(react$persistence_marked))) # selected homology classes idcs
      isolate({
        if ((length(idcs_selected) > 0 && !is.null(react$representations)) || react$dendrogram_redraw_zoom) {
          dendrogram_plotted <- TRUE # if TRUE and SVG export button was clicked, also export the SVG
          simplices_selected            <- react$persistence$inds$death[idcs_selected]
          react$representations.reduced <- lapply(react$representations, function(x) x[x %in% simplices_selected])
          
          if (((!is.null(react$dendrogram_marked_leaves.A) || !is.null(react$dendrogram_marked_leaves.B)) &&
               (react$dendrogram_redraw_highlights || react$image_export.dendrogram)) || react$dendrogram_redraw_zoom) {
            dendrogram_plot               <- trajectories_dendrogram(precomputed_dendrogram        = react$dendrogram,
                                                                     precomputed_dendrogram_labels = react$dendrogram_labels,
                                                                     leaves_to_highlight.A         = highlights.A,
                                                                     leaves_to_highlight.B         = highlights.B,
                                                                     leaves_to_highlight.zoom      = dendrogram_zoom_idcs)
            react$dendrogram_redraw_zoom       <- FALSE
            react$dendrogram_redraw_highlights <- FALSE
          } else {
            dendrogram_plot               <- trajectories_dendrogram(pers           = react$persistence,
                                                                     repre.reduced  = react$representations.reduced,
                                                                     perc           = perc,
                                                                     out.dendrogram = react$dendrogram,
                                                                     out.data       = react$dendrogram_data,
                                                                     out.classif    = react$dendrogram_classes,
                                                                     out.labels     = react$dendrogram_labels)
          }
          plot(dendrogram_plot)
        }
      })
      if (dendrogram_plotted && react$image_export.dendrogram) {
        if (react$image_export_format == 'SVG') {
          svg(filename = paste0('Dendrogram_', Sys.time(), '.svg'))
        } else {
          png(filename = paste0('Dendrogram_', Sys.time(), '.png'), width = 480, height = 5400)
        }
        plot(dendrogram_plot)
        dev.off()
        react$dendrogram_redraw_highlights <- TRUE
        react$image_export.dendrogram <- FALSE
      }
    }
    if (!dendrogram_plotted) {
      .draw_placeholder()
    }
  })
  
  output$plot_dendrogram_zoom <- renderPlot({
    session$resetBrush('selector_dendrogram')
    if (!react$dendrogram_zoom_active) {
      react$dendrogram_zoom_redraw_highlights <- TRUE
    }
    react$dendrogram_zoom_active <- TRUE
    par(mar = c(0, 0, 0, 0))
    perc <- react$dendrogram_leaf_perc_cutoff
    dendrogram_zoom_plotted <- FALSE
    dendrogram_zoom_plot    <- NULL
    highlights.A       <- react$dendrogram_marked_leaves.A
    highlights.B       <- react$dendrogram_marked_leaves.B
    zoom_idcs          <- react$dendrogram_zoom_idcs
    #if (react$dendrogram_zoom_ready) {
    if (react$dendrogram_draw_zoomed_in) {
      idcs_selected <- as.numeric(unlist(rownames(react$persistence_marked))) # selected homology classes idcs
      isolate({
        if ((length(idcs_selected) > 0 && !is.null(react$representations)) || react$dendrogram_redraw_zoom) {
          dendrogram_zoom_plotted       <- TRUE # if TRUE and SVG export button was clicked, also export the SVG
          simplices_selected            <- react$persistence$inds$death[idcs_selected]
          react$representations.reduced <- lapply(react$representations, function(x) x[x %in% simplices_selected])
          
          
          if (react$dendrogram_zoom_ready && (!is.null(react$dendrogram_marked_leaves.A) || !is.null(react$dendrogram_marked_leaves.B) &&
                                              (react$dendrogram_zoom_redraw_highlights || react$image_export.dendrogram_zoom))) {
            dendrogram_zoom_plot        <- trajectories_dendrogram(  precomputed_dendrogram        = react$dendrogram_zoom,
                                                                     precomputed_dendrogram_labels = react$dendrogram_zoom_labels,
                                                                     leaves_to_highlight.A         = highlights.A,
                                                                     leaves_to_highlight.B         = highlights.B,
                                                                     zoom_idcs                     = zoom_idcs)
            react$dendrogram_zoom_ready <- TRUE
            react$dendrogram_zoom_redraw_highlights <- FALSE
          } else {
            dendrogram_zoom_plot        <- trajectories_dendrogram(  pers           = react$persistence,
                                                                     repre.reduced  = react$representations.reduced,
                                                                     perc           = perc,
                                                                     zoom_idcs      = zoom_idcs,
                                                                     leaves_to_highlight.A         = highlights.A,
                                                                     leaves_to_highlight.B         = highlights.B,
                                                                     out.dendrogram = react$dendrogram_zoom,
                                                                     out.data       = react$dendrogram_zoom_data,
                                                                     out.classif    = react$dendrogram_zoom_classes,
                                                                     out.labels     = react$dendrogram_zoom_labels)
            react$dendrogram_zoom_ready <- TRUE
            if (!is.null(react$dendrogram_marked_leaves.A) || !is.null(react$dendrogram_marked_leaves.B)) {
              react$dendrogram_zoom_redraw_highlights <- TRUE
            }
          }
          plot(dendrogram_zoom_plot)
        }
      })
    }
    if (!dendrogram_zoom_plotted) {
      .draw_placeholder(picture = 'petal')
    }
  })
  
  # Leaf magnitude cutoff slider (in percentages)
  observeEvent(input$slider_dendrogram_leaf_cutoff, {
    react$dendrogram_zoom       <- NA
    react$dendrogram_zoom_ready <- FALSE
    react$dendrogram_zoom_idcs  <- NULL
    session$resetBrush('selector_dendrogram_zoom')
    react$dendrogram_draw_zoomed_in <- FALSE
    react$dendrogram_leaf_perc_cutoff <- input$slider_dendrogram_leaf_cutoff
  })
  
  # Buttons & logs
  observeEvent((input$btn_dendrogram_mark_leaves | input$btn_dendrogram_zoom_mark_leaves), {
    if (react$trajectories_group == 'A') {
      react$trajectories_marked.A        <- unique(unlist(c(react$trajectories_marked.A,        react$dendrogram_selection)))
      react$dendrogram_marked_leaves.A   <- unique(unlist(c(react$dendrogram_marked_leaves.A,   react$dendrogram_selected_leaves)))
      react$dendrogram_redraw_highlights <- TRUE
      
      react$trajectories_marked_idcs.A   <- unique(unlist(c(react$trajectories_marked_idcs.A),  react$dendrogram_selected_idcs))
    } else if (react$trajectories_group == 'B') {
      react$trajectories_marked.B        <- unique(unlist(c(react$trajectories_marked.B,        react$dendrogram_selection)))
      react$dendrogram_marked_leaves.B   <- unique(unlist(c(react$dendrogram_marked_leaves.B, react$dendrogram_selected_leaves)))
      react$dendrogram_redraw_highlights <- TRUE
      
      react$trajectories_marked_idcs.B   <- unique(unlist(c(react$trajectories_marked_idcs.B),  react$dendrogram_selected_idcs))
    }
  })
  
  observeEvent((input$btn_dendrogram_clear_marked_leaves | input$btn_dendrogram_zoom_clear_marked_leaves), {
    if (react$trajectories_group == 'A') {
      if (!is.null(react$trajectories_marked.A) && is.null(react$trajectories_marked.B) && !is.null(react$dendrogram_zoom_idcs)) {
        react$dendrogram_redraw_zoom <- TRUE
      }
      react$trajectories_marked.A          <- react$dendrogram_marked_leaves.A <- react$trajectories_marked_idcs.A <- NULL
      react$trajectories_junk.A            <- integer(0)
      react$tracked_markers_last_removed.A <- NULL
      react$dendrogram_redraw_highlights   <- TRUE
    } else if (react$trajectories_group == 'B') {
      if (!is.null(react$trajectories_marked.B) && is.null(react$trajectories_marked.A) && !is.null(react$dendrogram_zoom_idcs)) {
        react$dendrogram_redraw_zoom <- TRUE
      }
      react$trajectories_marked.B          <- react$dendrogram_marked_leaves.B <- react$trajectories_marked_idcs.B <- NULL
      react$trajectories_junk.B            <- integer(0)
      react$tracked_markers_last_removed.B <- NULL
      react$dendrogram_redraw_highlights   <- TRUE
    }
  })
  
  pin_batch_modal <- function(default_name, failed = FALSE) {
    modalDialog(
      textInput('input_pin_batch_name', 'Pinned batch name', placeholder = '', value = default_name),
      if (failed) div(tags$b('Field is empty', style = 'colour: red')),
      footer = tagList(
        modalButton('Cancel'),
        actionButton('btn_pin_batch_name', 'Confirm')
      )
    )
  }
  
  observeEvent(react$trajectories_to_pin, {
    if (is.null(react$output_ff)) {
      if (is.null(event_sel)) {
        layout_X <- layout.df[[react$layout_name]]$X * 100
        layout_Y <- layout.df[[react$layout_name]]$Y * 100
        ps       <- as.numeric(as.factor(react$pseudotime$res))
      } else {
        layout_X <- layout_Y <- rep(-100, nrow(input_ff))
        layout_X[event_sel] <- layout.df[[react$layout_name]]$X * 100
        layout_Y[event_sel] <- layout.df[[react$layout_name]]$Y * 100
        ps           <- rep(-100, nrow(input_ff))
        ps[event_sel] <- as.numeric(as.factor(react$pseudotime$res))
      }
      react$output_ff <- fcs.add_col(
        fcs.add_col(
          fcs.add_col(
            input_ff, layout_X, colname = 'dimension_reduction_1'
          ),
          layout_Y, colname = 'dimension_reduction_2'
        ),
        ps,
        colname = 'total_pseudotime'
      )
    }
    
    react$trajectories_pinned_batches_count <- react$trajectories_pinned_batches_count + 1
    
    showModal(pin_batch_modal(default_name = paste0('pathway_batch_', sprintf("%03d", react$trajectories_pinned_batches_count))))
  })
  
  observeEvent(input$btn_pin_batch_name, {
    
    if (input$input_pin_batch_name != '') {
      walks        <- list()
      walks$v      <- unlist(react$trajectories_random_walks)
      walks$starts <- c(1, cumsum(sapply(react$trajectories_random_walks, length)) + 1)
      react$output_ff <- .add_path_info_to_fcs(ff                       = react$output_ff,
                                               pseudotime               = react$pseudotime,
                                               all_walks                = walks,
                                               trajectories_of_interest = react$trajectories_to_pin,
                                               id                       = input$input_pin_batch_name,
                                               event_sel                = event_sel)
      removeModal()
    } else {
      showModal(pin_batch_modal(failed = TRUE))
    }
  })
  
  observeEvent(react$trajectories_pinned_batches_count, {
    count <- react$trajectories_pinned_batches_count
    output$log_pinned_batches_count <- output$pinned_batches_count_zoom <- renderPrint({
      cat(paste(count, ' ', if (count == 1) { 'batch' } else { 'batches' }, 'pinned'), '\n')
    })
  })
  
  export_fcs_modal <- function(failed = FALSE) {
    modalDialog(
      textInput('input_export_fcs_name', 'Output FCS file name', placeholder = '', value = 'output.FCS'),
      if (failed) div(tags$b('Field is empty', style = 'colour: red')),
      footer = tagList(
        modalButton('Cancel'),
        actionButton('btn_export_fcs_modal_save', 'Save')
      )
    )
  }
  
  
  
  observeEvent((input$btn_trajectories_export_fcs), {
    if (is.null(react$output_ff)) {
      if (is.null(event_sel)) {
        layout_X <- layout.df[[react$layout_name]]$X * 100
        layout_Y <- layout.df[[react$layout_name]]$Y * 100
        ps       <- as.numeric(as.factor(react$pseudotime$res))
      } else {
        layout_X <- layout_Y <- rep(-100, nrow(input_ff))
        layout_X[event_sel] <- layout.df[[react$layout_name]]$X * 100
        layout_Y[event_sel] <- layout.df[[react$layout_name]]$Y * 100
        ps           <- rep(-100, nrow(input_ff))
        ps[event_sel] <- as.numeric(as.factor(react$pseudotime$res))
      }
      react$output_ff <- fcs.add_col(
        fcs.add_col(
          fcs.add_col(
            input_ff, layout_X, colname = 'dimension_reduction_1'
          ),
          layout_Y, colname = 'dimension_reduction_2'
        ),
        ps,
        colname = 'total_pseudotime'
      )
    }
    showModal(export_fcs_modal())
  })
  
  observeEvent((input$btn_trajectories_clear_pinned_trajectories), {
    if (is.null(react$output_ff)) {
      if (is.null(event_sel)) {
        layout_X <- layout.df[[react$layout_name]]$X * 100
        layout_Y <- layout.df[[react$layout_name]]$Y * 100
        ps       <- react$pseudotime$res
      } else {
        layout_X <- layout_Y <- rep(-100, nrow(input_ff))
        layout_X[event_sel] <- layout.df[[react$layout_name]]$X * 100
        layout_Y[event_sel] <- layout.df[[react$layout_name]]$Y * 100
        ps           <- rep(-100, nrow(input_ff))
        ps[event_sel] <- react$pseudotime$res
      }
      react$output_ff <- fcs.add_col(
        fcs.add_col(
          fcs.add_col(
            input_ff, layout_X, colname = 'dimension_reduction_1'
          ),
          layout_Y, colname = 'dimension_reduction_2'
        ),
        ps,
        colname = 'total_pseudotime'
      )
    }
    react$trajectories_pinned_batches_count <- 0
    react$trajectories_pinned               <- NULL
  })
  observeEvent(input$btn_export_fcs_modal_save, {
    if (input$input_export_fcs_name != '') {
      flowCore::write.FCS(react$output_ff, input$input_export_fcs_name)
      removeModal()
    } else {
      showModal(export_fcs_modal(failed = TRUE))
    }
  })
  
  observeEvent((input$btn_dendrogram_export_image), {
    react$image_export.dendrogram <- TRUE
  })
  
  observeEvent(input$btn_dendrogram_zoom, {
    df  <- data.frame(y = seq(0, 1, length.out = length(react$dendrogram_classes)))
    pts <- brushedPoints(df, input$selector_dendrogram, yvar = 'y')
    pts <- rev(as.numeric(rownames(pts)))
    if (length(pts) < 2) {
      react$dendrogram_zoom_idcs      <- NULL
      react$dendrogram_zoom_ready     <- FALSE
      react$dendrogram_redraw_zoom    <- TRUE
      react$dendrogram_draw_zoomed_in <- FALSE
    } else {
      X      <- react$dendrogram_data$labels[pts, 1]
      react$dendrogram_zoom_idcs      <- c(min(X), max(X))
      react$dendrogram_zoom_ready     <- FALSE
      react$dendrogram_redraw_zoom    <- TRUE
      react$dendrogram_draw_zoomed_in <- TRUE
    }
  })
  
  # observeEvent((input$btn_dendrogram_mark_leaves | input$btn_dendrogram_zoom_mark_leaves), {
  #   if (react$trajectories_group == 'A') {
  #     react$trajectories_marked.A        <- unique(unlist(c(react$trajectories_marked.A,        react$dendrogram_selection)))
  #     react$dendrogram_marked_leaves.A   <- unique(unlist(c(react$dendrogram_marked_leaves.A,   react$dendrogram_selected_leaves)))
  #     react$dendrogram_redraw_highlights <- TRUE
  #     react$trajectories_marked_idcs.A   <- unique(unlist(c(react$trajectories_marked_idcs.A),  react$dendrogram_selected_idcs))
  #   } else if (react$trajectories_group == 'B') {
  #     react$trajectories_marked.B        <- unique(unlist(c(react$trajectories_marked.B,        react$dendrogram_selection)))
  #     react$dendrogram_marked_leaves.B   <- unique(unlist(c(react$dendrogram_marked_leaves.B, react$dendrogram_selected_leaves)))
  #     react$dendrogram_redraw_highlights <- TRUE
  #     react$trajectories_marked_idcs.B   <- unique(unlist(c(react$trajectories_marked_idcs.B),  react$dendrogram_selected_idcs))
  #   }
  # })
  
  
  observeEvent(input$btn_trajectories_group, {
    react$trajectories_group <- input$btn_trajectories_group
  })
  
  observeEvent(input$btn_trajectories_pin_trajectories.A, {
    react$trajectories_to_pin <- sort(unique(react$trajectories_marked.A))
  })
  
  observeEvent(input$btn_trajectories_pin_trajectories.B, {
    react$trajectories_to_pin <- sort(unique(react$trajectories_marked.B))
  })
  
  output$log_dendrogram_selected <- renderPrint({
    if (any(c(!is.na(react$dendrogram), !is.na(react$dendrogram_zoom)))) {
      if (react$dendrogram_zoom_active) {
        df      <- data.frame(y = seq(0, 1, length.out = length(react$dendrogram_zoom_classes)))
        pts     <- brushedPoints(df, input$selector_dendrogram_zoom, yvar = 'y')
        pts     <- rev(as.numeric(rownames(pts)))
        leaves  <- names(react$dendrogram_zoom_classes)[pts]
        #uniques                          <- which(!duplicated(leaves))
        react$dendrogram_selected_leaves <- leaves#[uniques]
        react$dendrogram_selected_idcs   <- as.vector(unlist(react$dendrogram_zoom_classes[pts]))#[uniques]
        react$dendrogram_selection       <- as.numeric(unlist(react$dendrogram_zoom_classes[pts]))#[uniques]
        sizes                            <- unlist(sapply(react$dendrogram_zoom_classes[pts], length))#[uniques]
      } else {
        df      <- data.frame(y = seq(0, 1, length.out = length(react$dendrogram_classes)))
        pts     <- brushedPoints(df, input$selector_dendrogram, yvar = 'y')
        pts     <- rev(as.numeric(rownames(pts)))
        leaves  <- names(react$dendrogram_classes)[pts]
        #uniques                          <- which(!duplicated(leaves))
        react$dendrogram_selected_leaves <- leaves#[uniques]
        react$dendrogram_selected_idcs   <- as.vector(unlist(react$dendrogram_classes[pts]))#[uniques]
        react$dendrogram_selection       <- as.numeric(unlist(react$dendrogram_classes[pts]))#[uniques]
        sizes                            <- unlist(sapply(react$dendrogram_classes[pts], length))#[uniques]
      }
      cat(sizes, sep = ', ')
    }
  })
  
  output$log_dendrogam_marked.A <- output$log_dendrogam_marked.B <- renderPrint({ cat('(none)') })
  
  observeEvent(react$trajectories_marked.A, {
    react$trajectories_junk.A <- react$trajectories_junk.A[!react$trajectories_junk.A %in% react$trajectories_marked.A]
    junk_count                <- length(react$trajectories_junk.A)
    marked_count              <- length(react$trajectories_marked.A)
    msg.erased                <- if (junk_count > 0) { paste0('(', junk_count, ' zapped)') } else { '' }
    output$log_trajectories_marked_count.A <- renderPrint({ cat('N =', marked_count, msg.erased) })
    if (marked_count > 250) {
      output$log_trajectories_marked.A <- renderPrint({ cat(paste(react$trajectories_marked.A[1:250], sep = ', '), paste('...and', marked_count - 250, 'more'), sep = ', ') })
    } else if (marked_count > 0) {
      output$log_trajectories_marked.A <- renderPrint({ cat(react$trajectories_marked.A, sep = ', ') })
    } else {
      output$log_trajectories_marked.A <- renderPrint({ cat('(none)\n') })
    }
  })
  
  observeEvent(react$trajectories_marked.B, {
    react$trajectories_junk.B <- react$trajectories_junk.B[!react$trajectories_junk.B %in% react$trajectories_marked.B]
    junk_count                <- length(react$trajectories_junk.B)
    marked_count              <- length(react$trajectories_marked.B)
    msg.erased                <- if (junk_count > 0) { paste0('(', junk_count, ' zapped)') } else { '' }
    output$log_trajectories_marked_count.B <- renderPrint({ cat('N =', marked_count, msg.erased) })
    if (marked_count > 250) {
      output$log_trajectories_marked.B <- renderPrint({ cat(paste(react$trajectories_marked.B[1:250], sep = ', '), paste('...and', marked_count - 250, 'more'), sep = ', ') })
    } else if (marked_count > 0) {
      output$log_trajectories_marked.B <- renderPrint({ cat(react$trajectories_marked.B, sep = ', ') })
    } else {
      output$log_trajectories_marked.B <- renderPrint({ cat('(none)\n') })
    }
  })
  
  ## Trajectories layout
  # Plot
  output$plot_layout_trajectories <- renderPlot({
    if (react$image_export.trajectories) {
      if (react$image_export_format == 'SVG') {
        svg(filename = paste0('Trajectories_', Sys.time(), '.svg'))
      } else {
        png(filename = paste0('Trajectories_', Sys.time(), '.png'), width = 1000, height = 900)
      }
      .plot_trajectories_full(layout[[react$layout_name]],
                              react$trajectories_random_walks,
                              react$trajectories_marked.A,
                              react$trajectories_marked.B,
                              flip_colours = react$layout_trajectories_flip_colours,
                              pseudotime_highlight_bounds = react$pseudotime_highlight_bounds,
                              pseudotime = if (is.null(react$pseudotime_highlight_bounds)) { NULL } else { react$pseudotime },
                              highlight_in_background = react$layout_trajectories_highlight_in_background,
                              selected_trajectory_points = react$selected_trajectory_points)
      dev.off()
      react$image_export.trajectories <- FALSE
    } else {
      par(mar = c(0, 0, 0, 0))
      if (!is.null(react$layout_name) && !is.null(layout[[react$layout_name]])) {
        .plot_trajectories(layout[[react$layout_name]],
                           react$trajectories_random_walks,
                           react$trajectories_marked.A,
                           react$trajectories_marked.B,
                           flip_colours = react$layout_trajectories_flip_colours,
                           pseudotime_highlight_bounds = react$pseudotime_highlight_bounds,
                           pseudotime = if (is.null(react$pseudotime_highlight_bounds)) { NULL } else { react$pseudotime },
                           highlight_in_background = react$layout_trajectories_highlight_in_background, selected_trajectory_points = react$selected_trajectory_points,
                           pointsize = 1)
      }
      
    }
  })
  
  # Buttons
  observeEvent(input$btn_layout_trajectories_remove_highlight, {
    react$pseudotime_highlight_bounds  <- NULL
  })
  observeEvent(input$btn_layout_trajectories_highlight_in_background, {
    react$layout_trajectories_highlight_in_background <- !react$layout_trajectories_highlight_in_background
  })
  observeEvent(input$btn_layout_trajectories_flip_colours, {
    react$layout_trajectories_flip_colours <- !react$layout_trajectories_flip_colours
  })
  observeEvent(input$btn_layout_trajectories_export_image, {
    react$image_export.trajectories <- TRUE
  })
  
  ## Marker expression trackers
  # Inputs
  observe({
    updateSelectInput(session, 'input_tracked_markers.A', choices = markers)
    updateSelectInput(session, 'input_tracked_markers.B', choices = markers)
  })
  observeEvent(input$input_tracked_markers.A, {
    react$tracked_markers.A <- input$input_tracked_markers.A
    if (is.null(react$tracked_markers.B)) {
      updateSelectInput(session, 'input_tracked_markers.B', selected = input$input_tracked_markers.A)
    }
  })
  observeEvent(input$input_tracked_markers.B, {
    react$tracked_markers.B <- input$input_tracked_markers.B
    if (is.null(react$tracked_markers.A)) {
      updateSelectInput(session, 'input_tracked_markers.A', selected = input$input_tracked_markers.B)
    }
  })
  observeEvent(input$input_trackers_scaling_exponent, {
    if (is.numeric(input$input_trackers_scaling_exponent)) {
      if (input$input_trackers_scaling_exponent > 0.001 && input$input_trackers_scaling_exponent < 5) {
        react$trackers_scaling_exponent <- input$input_trackers_scaling_exponent
      }
    }
  })
  observeEvent(input$input_trackers_n_segments, {
    if (is.numeric(input$input_trackers_n_segments)) {
      if (input$input_trackers_n_segments > 3 && input$input_trackers_n_segments < 1000) {
        react$trackers_n_segments <- as.integer(input$input_trackers_n_segments)
      }
    }
  })
  observeEvent(input$check_trackers_large_base_size, {
    react$trackers_large_base_size <- input$check_trackers_large_base_size
  })
  observeEvent(input$btn_tracked_markers_remove_trajectories.A, {
    if (length(react$tracked_markers.A) == 1) {
      pts <- brushedPoints(react$tracked_markers_stats.A,
                           input$selector_tracked_markers.A,
                           xvar = 'segment', yvar = 'expression')
      n_selected <- length(unique(pts$walk))
      n_total <- length(unlist(react$trajectories_marked.A))
      if (n_selected < n_total) {
        junk                                 <- react$trajectories_marked.A[unique(pts$walk)]
        react$tracked_markers_last_removed.A <- junk
        react$trajectories_junk.A            <- na.omit(unique(c(react$trajectories_junk.A, junk)))
        nonjunk                              <- !react$trajectories_marked.A %in% junk
        react$trajectories_marked.A          <- na.omit(react$trajectories_marked.A[nonjunk])
        session$resetBrush('selector_tracked_markers.A')
      }
    }
  })
  observeEvent(input$btn_tracked_markers_remove_trajectories.B, {
    if (length(react$tracked_markers.B) == 1) {
      pts <- brushedPoints(react$tracked_markers_stats.B,
                           input$selector_tracked_markers.B,
                           xvar = 'segment', yvar = 'expression')
      n_selected <- length(unique(pts$walk))
      n_total <- length(unlist(react$trajectories_marked.B))
      if (n_selected < n_total) {
        junk                                 <- react$trajectories_marked.B[unique(pts$walk)]
        react$tracked_markers_last_removed.B <- junk
        react$trajectories_junk.B            <- na.omit(unique(c(react$trajectories_junk.B, junk)))
        nonjunk                              <- !react$trajectories_marked.B %in% junk
        react$trajectories_marked.B          <- na.omit(react$trajectories_marked.B[nonjunk])
        session$resetBrush('selector_tracked_markers.B')
      }
    }
  })
  
  observeEvent(input$btn_tracked_markers_undo_remove_trajectories.A, {
    if (!is.null(react$tracked_markers_last_removed.A)) {
      react$trajectories_marked.A          <- unique(unlist(c(react$trajectories_marked.A, react$tracked_markers_last_removed.A)))
      idcs                                 <- which(react$trajectories_junk.A %in% react$tracked_markers_last_removed.A)
      react$trajectories_junk.A            <- react$trajectories_junk.A[-idcs]
      react$tracked_markers_last_removed.A <- NULL
    }
  })
  observeEvent(input$btn_tracked_markers_undo_remove_trajectories.B, {
    if (!is.null(react$tracked_markers_last_removed.B)) {
      react$trajectories_marked.B          <- unique(unlist(c(react$trajectories_marked.B, react$tracked_markers_last_removed.B)))
      idcs                                 <- which(react$trajectories_junk.B %in% react$tracked_markers_last_removed.B)
      react$trajectories_junk.B            <- react$trajectories_junk.B[-idcs]
      react$tracked_markers_last_removed.B <- NULL
    }
  })
  
  observeEvent(input$btn_tracked_markers_highlight_segments.A, {
    if (!is.null(react$tracked_markers.A)) {
      pts <- brushedPoints(react$tracked_markers_stats.A,
                           input$selector_tracked_markers.A,
                           xvar = 'segment', yvar = 'expression')
      
      if (nrow(pts)>0){
        
        ##----MODIFIED by JS - temp
        selected_trajectory_points<-NULL
        pts$inds_char<-as.character(pts$inds_char)
        for (i in 1:nrow(pts)){
          
          if (pts$inds_char[i]!="NULL") selected_trajectory_points<-c(selected_trajectory_points,as.numeric(unlist(strsplit(pts$inds_char[i],split=","))))
        }
        react$selected_trajectory_points<-unique(selected_trajectory_points)
      }
      
      ##--END of JS
      
      
      if (nrow(pts) > 0) {
        highlighted_segments <- sort(unique(pts$segment))
        react$pseudotime_highlight_bounds <- react$tracked_markers_pseudotime_bounds.A[c(min(highlighted_segments), max(highlighted_segments) + 1)]
      } else {
        react$pseudotime_highlight_bounds <- NULL
      }
      session$resetBrush('selector_tracked_markers.A')
    }
  })
  
  observeEvent(input$btn_tracked_markers_highlight_segments.B, {
    if (!is.null(react$tracked_markers.B)) {
      pts <- brushedPoints(react$tracked_markers_stats.B,
                           input$selector_tracked_markers.B,
                           xvar = 'segment', yvar = 'expression')
      
      if (nrow(pts)>0){
        
        ##----MODIFIED by JS - temp
        selected_trajectory_points<-NULL
        pts$inds_char<-as.character(pts$inds_char)
        for (i in 1:nrow(pts)){
          
          if (pts$inds_char[i]!="NULL") selected_trajectory_points<-c(selected_trajectory_points,as.numeric(unlist(strsplit(pts$inds_char[i],split=","))))
        }
        react$selected_trajectory_points<-unique(selected_trajectory_points)
      }
      
      ##--END of JS
      
      if (nrow(pts) > 0) {
        highlighted_segments <- sort(unique(pts$segment))
        react$pseudotime_highlight_bounds <- react$tracked_markers_pseudotime_bounds.B[c(min(highlighted_segments), max(highlighted_segments) + 1)]
      } else {
        react$pseudotime_highlight_bounds <- NULL
      }
      session$resetBrush('selector_tracked_markers.B')
    }
  })
  
  observeEvent(input$btn_tracked_populations_highlight_segments.A, {
    if (!is.null(react$tracked_populations.A)) {
      pts <- brushedPoints(react$tracked_populations_stats.A,
                           input$selector_tracked_populations.A,
                           xvar = 'segment', yvar = 'count')
      if (nrow(pts)>0){
        selected_trajectory_points<-NULL
        pts$inds_char<-as.character(pts$inds_char)
        for (i in 1:nrow(pts)){
          
          if (pts$inds_char[i]!="NULL") selected_trajectory_points<-c(selected_trajectory_points,as.numeric(unlist(strsplit(pts$inds_char[i],split=","))))
        }
        react$selected_trajectory_points<-unique(selected_trajectory_points)
        
        highlighted_segments <- sort(unique(pts$segment))
        react$pseudotime_highlight_bounds <- react$tracked_markers_pseudotime_bounds.B[c(min(highlighted_segments), max(highlighted_segments) + 1)]
      } else {
        react$pseudotime_highlight_bounds <- NULL
      }
      session$resetBrush('selector_tracked_populations.A')
    }
  })
  
  observeEvent(input$btn_tracked_populations_highlight_segments.B, {
    if (!is.null(react$tracked_populations.B)) {
      pts <- brushedPoints(react$tracked_populations_stats.B,
                           input$selector_tracked_populations.B,
                           xvar = 'segment', yvar = 'count')
      if (nrow(pts)>0){
        
        ##----MODIFIED by JS - temp
        selected_trajectory_points<-NULL
        pts$inds_char<-as.character(pts$inds_char)
        for (i in 1:nrow(pts)){
          
          if (pts$inds_char[i]!="NULL") selected_trajectory_points<-c(selected_trajectory_points,as.numeric(unlist(strsplit(pts$inds_char[i],split=","))))
        }
        react$selected_trajectory_points<-unique(selected_trajectory_points)
        
        highlighted_segments <- sort(unique(pts$segment))
        react$pseudotime_highlight_bounds <- react$tracked_markers_pseudotime_bounds.B[c(min(highlighted_segments), max(highlighted_segments) + 1)]
      } else {
        react$pseudotime_highlight_bounds <- NULL
      }
      session$resetBrush('selector_tracked_populations.B')
    }
  })
  
  observeEvent(input$btn_tracked_markers_export_image.A, {
    react$image_export.tracked_markers.A <- TRUE
  })
  observeEvent(input$btn_tracked_markers_export_image.B, {
    react$image_export.tracked_markers.B <- TRUE
  })
  observeEvent(input$btn_tracked_populations_export_image.A, {
    react$image_export.tracked_populations.A <- TRUE
  })
  observeEvent(input$btn_tracked_populations_export_image.B, {
    react$image_export.tracked_populations.B <- TRUE
  })
  
  # Plots
  output$plot_tracked_markers.A <- renderPlot({
    session$resetBrush('selector_tracked_markers.A')
    if (!is.null(react$trajectories_marked.A) && !is.null(react$tracked_markers.A)) {
      react$tracked_markers_ready.A <- TRUE
      p <- .plot_tracked_markers(react$trajectories_random_walks,
                                 react$trajectories_marked.A,
                                 tv,
                                 react$pseudotime,
                                 markers  = react$tracked_markers.A,
                                 n.part   = react$trackers_n_segments,
                                 exp.part = react$trackers_scaling_exponent,
                                 large_base_size = react$trackers_large_base_size,
                                 grey = !react$image_export.tracked_markers.A)
      react$tracked_markers_stats.A <- p$stats
      react$tracked_markers_pseudotime_bounds.A <- p$pseudotime_bounds
      if (react$image_export.tracked_markers.A) {
        if (react$image_export_format == 'SVG') {
          svg(filename = paste0('Markers_A_', Sys.time(), '.svg'), width = 10, height = 6)
        } else {
          png(filename = paste0('Markers_A_', Sys.time(), '.png'), width = 750, height = 480)
        }
        pp <- p$plot
        if (react$image_export_format != 'SVG') {
          pp <- pp +
            theme(text          = element_text(size = 24),
                  plot.title    = element_text(size = 20),
                  plot.subtitle = element_text(size = 19),
                  legend.title  = element_text(size = 24),
                  legend.text   = element_text(size = 22)) +
            geom_line(size = 2)
        }
        plot(pp)
        dev.off()
        react$image_export.tracked_markers.A <- FALSE
      } else {
        p$plot
      }
    } else {
      react$tracked_markers_ready.A <- FALSE
    }
  })
  output$plot_tracked_markers.B <- renderPlot({
    session$resetBrush('selector_tracked_markers.B')
    if (!is.null(react$trajectories_marked.B) && !is.null(react$tracked_markers.B)) {
      react$tracked_markers_ready.B <- TRUE
      p <- .plot_tracked_markers(react$trajectories_random_walks,
                                 react$trajectories_marked.B,
                                 tv,
                                 react$pseudotime,
                                 markers  = react$tracked_markers.B,
                                 n.part   = react$trackers_n_segments,
                                 exp.part = react$trackers_scaling_exponent,
                                 large_base_size = react$trackers_large_base_size,
                                 grey = !react$image_export.tracked_markers.B)
      react$tracked_markers_stats.B <- p$stats
      react$tracked_markers_pseudotime_bounds.B <- p$pseudotime_bounds
      if (react$image_export.tracked_markers.B) {
        if (react$image_export_format == 'SVG') {
          svg(filename = paste0('Markers_B_', Sys.time(), '.svg'), width = 10, height = 6)
        } else {
          png(filename = paste0('Markers_B_', Sys.time(), '.png'), width = 750, height = 480)
        }
        pp <- p$plot
        if (react$image_export_format != 'SVG') {
          pp <- pp +
            theme(text          = element_text(size = 24),
                  plot.title    = element_text(size = 20),
                  plot.subtitle = element_text(size = 19),
                  legend.title  = element_text(size = 24),
                  legend.text   = element_text(size = 22)) +
            geom_line(size = 2)
        }
        plot(pp)
        dev.off()
        react$image_export.tracked_markers.B <- FALSE
      } else {
        p$plot
      }
    } else {
      react$tracked_markers_ready.B <- FALSE
    }
  })
  
  ## Population composition trackers
  # Inputs
  observeEvent(input$check_tracked_populations_log2_transform, {
    react$tracked_populations_log2_transform <- input$check_tracked_populations_log2_transform
  })
  observe({
    updateSelectInput(session, 'input_tracked_populations.A', choices = labels.aligned[[react$labels_name]])
    updateSelectInput(session, 'input_tracked_populations.B', choices = labels.aligned[[react$labels_name]])
  })
  observeEvent(input$input_tracked_populations.A, {
    react$tracked_populations.A <- input$input_tracked_populations.A
    if (is.null(react$tracked_populations.B)) {
      updateSelectInput(session, 'input_tracked_populations.B', selected = input$input_tracked_populations.A)
    }
  })
  observeEvent(input$input_tracked_populations.B, {
    react$tracked_populations.B <- input$input_tracked_populations.B
    if (is.null(react$tracked_populations.A)) {
      updateSelectInput(session, 'input_tracked_populations.A', selected = input$input_tracked_populations.B)
    }
  })
  
  # Plots
  output$plot_tracked_populations.A <- renderPlot({
    session$resetBrush('selector_tracked_populations.A')
    if (!is.null(react$trajectories_marked.A) && !is.null(react$tracked_populations.A)) {
      react$tracked_markers_ready.A <- TRUE
      p <- .plot_tracked_populations(react$trajectories_random_walks,
                                     react$trajectories_marked.A,
                                     tv,
                                     labels_name = react$labels_name,
                                     react$pseudotime,
                                     populations    = react$tracked_populations.A,
                                     n.part         = react$trackers_n_segments,
                                     exp.part       = react$trackers_scaling_exponent,
                                     log2_transform = react$tracked_populations_log2_transform,
                                     large_base_size = react$trackers_large_base_size,
                                     grey = !react$image_export.tracked_populations.A)
      react$tracked_populations_stats.A         <- p$stats
      react$tracked_markers_pseudotime_bounds.A <- p$pseudotime_bounds
      
      if (react$image_export.tracked_populations.A) {
        if (react$image_export_format == 'SVG') {
          svg(filename = paste0('Populations_A_', Sys.time(), '.svg'), width = 10, height = 6)
        } else {
          png(filename = paste0('Populations_A_', Sys.time(), '.png'), width = 750, height = 480)
        }
        pp <- p$plot
        if (react$image_export_format != 'SVG') {
          pp <- pp +
            theme(text          = element_text(size = 24),
                  plot.title    = element_text(size = 20),
                  plot.subtitle = element_text(size = 19),
                  legend.title  = element_text(size = 24),
                  legend.text   = element_text(size = 22)) +
            geom_line(size = 2)
        }
        plot(pp)
        dev.off()
        react$image_export.tracked_populations.A <- FALSE
      } else {
        p$plot
      }
    } else {
      react$tracked_markers_ready.A <- FALSE
    }
  })
  output$plot_tracked_populations.B <- renderPlot({
    session$resetBrush('selector_tracked_populations.B')
    if (!is.null(react$trajectories_marked.B) && !is.null(react$tracked_populations.B)) {
      react$tracked_markers_ready.B <- TRUE
      p <- .plot_tracked_populations(react$trajectories_random_walks,
                                     react$trajectories_marked.B,
                                     tv,
                                     labels_name = react$labels_name,
                                     react$pseudotime,
                                     populations    = react$tracked_populations.B,
                                     n.part         = react$trackers_n_segments,
                                     exp.part       = react$trackers_scaling_exponent,
                                     log2_transform = react$tracked_populations_log2_transform,
                                     large_base_size = react$trackers_large_base_size,
                                     grey = !react$image_export.tracked_populations.B)
      react$tracked_populations_stats.B         <- p$stats
      react$tracked_markers_pseudotime_bounds.B <- p$pseudotime_bounds
      
      if (react$image_export.tracked_populations.B) {
        if (react$image_export_format == 'SVG') {
          svg(filename = paste0('Populations_B_', Sys.time(), '.svg'), width = 10, height = 6)
        } else {
          png(filename = paste0('Populations_B_', Sys.time(), '.png'), width = 750, height = 480)
        }
        pp <- p$plot
        if (react$image_export_format != 'SVG') {
          pp <- pp +
            theme(text          = element_text(size = 24),
                  plot.title    = element_text(size = 20),
                  plot.subtitle = element_text(size = 19),
                  legend.title  = element_text(size = 24),
                  legend.text   = element_text(size = 22)) +
            geom_line(size = 2)
        }
        plot(pp)
        dev.off()
        react$image_export.tracked_populations.B <- FALSE
      } else {
        p$plot
      }
    } else {
      react$tracked_markers_ready.B <- FALSE
    }
  })
  ## TOOLTIPS ON HOVER
  addTooltip(session, 
             id    = 'btn_termini_mark_termini',
             title = 'Mark selected terminal nodes',
             trigger = 'hover',
             options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_termini_clear_termini',
             title = 'Clear marked terminal nodes',
             trigger = 'hover',
             options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_termini_update_walks_by_termini',
             title = 'Classify walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_termini_export_image',
             title = 'Export image of the layout above',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_left_show_gating',
             title = 'Show layout with labelled populations',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_termini_show_gating',
             title = 'Show layout with labelled populations',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_persistence_mark_classes',
             title = 'Mark selected homology classes',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_persistence_clear_classes',
             title = 'Clear marked homology classes',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_persistence_export_image',
             title = 'Export image of the above persistence diagram',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'switch_persistence_ratio',
             title = 'Show ratio of Death and Birth filtration values on y-axis, for data containing regions with dissimilar densities of points',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_dendrogram_mark_leaves',
             title = 'Mark selected leaves of HCP dendrogram',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_dendrogram_clear_marked_leaves',
             title = 'Clear marked leaves of the HCP dendrogram',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_dendrogram_zoom',
             title = 'Generate a zoomed-in view of the selected area',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_trajectories_export_fcs',
             title = 'Export FCS file with 2-d layout and all pinned walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_trajectories_clear_pinned_trajectories',
             title = 'Clear pinned walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_dendrogram_export_image',
             title = 'Export image of the HCP dendrogram',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_dendrogram_zoom_mark_leaves',
             title = 'Mark selected leaves of HCP dendrogram',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_dendrogram_zoom_mark_leaves',
             title = 'Mark selected leaves of HCP dendrogram',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'slider_dendrogram_leaf_cutoff',
             title = 'Set a threshold minimum percentage of all simulated walks per branch for pruning the HCP dendrogram',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_trajectories_pin_trajectories.A',
             title = 'Pin walks in group A',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_trajectories_pin_trajectories.B',
             title = 'Pin walks in group A',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_layout_trajectories_remove_highlight',
             title = 'Remove pseudotime segment highlighting in the above layout',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_layout_trajectories_highlight_in_background',
             title = 'Toggle pseudotime segment highlighting in background/foreground with respect to plotted walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_layout_trajectories_flip_colours',
             title = 'Toggle plotting group A walks in foreground/background with respect to group B walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_layout_trajectories_export_image',
             title = 'Export image of the above layout',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'input_trackers_scaling_exponent',
             title = 'Scaling exponent for pseudotime segmentation',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'input_trackers_n_segments',
             title = 'Number of pseudotime segments',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'check_trackers_large_base_size',
             title = 'Use larger font for in the plots below',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_remove_trajectories.A',
             title = 'Remove all simulated walks passing through selection from group A',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_remove_trajectories.B',
             title = 'Remove all simulated walks passing through selection from group B',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_undo_remove_trajectories.A',
             title = 'Undo last removal walks in group A',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_undo_remove_trajectories.B',
             title = 'Undo last removal walks in group B',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_highlight_segments.A',
             title = 'Highlight selected pseudotime segments in the layout above',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_highlight_segments.B',
             title = 'Highlight selected pseudotime segments in the layout above',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_export_image.A',
             title = 'Export image of tracked marker expression changes in group A of simulated walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_markers_export_image.B',
             title = 'Export image of tracked marker expression changes in group B of simulated walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_populations_highlight_segments.A',
             title = 'Highlight selected pseudotime segments in the layout above',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_populations_highlight_segments.B',
             title = 'Highlight selected pseudotime segments in the layout above',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_populations_export_image.A',
             title = 'Export image of tracked representation of populations along pseudotime progression in group A of simulated walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_tracked_populations_export_image.B',
             title = 'Export image of tracked representation of populations along pseudotime progression in group B of simulated walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'check_tracked_populations_log2_transform',
             title = 'Apply log2 transform to sizes of labelled populations for tracking representation of populations along pseudotime progression (good for extremely abundant and rare populations)',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_trackers_enlarge',
             title = 'Display larger marker and population trackers for the two groups of simulated walks',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
  addTooltip(session, 
             id    = 'btn_help',
             title = 'Show instructions for new users',
             trigger = 'hover',     options = list(delay = list(show=500))
  )
}