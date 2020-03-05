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
  layout           <- tv$layout
  event_sel        <- tv$events_sel # selected events from input FCS file
  markers          <- colnames(tv$data)
  labels.unique <- unique(tv$labels)
  
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
  # Terminal nodes selection
  react$termini_selected             <- NULL # selected terminal nodes idcs
  react$termini_marked               <- NULL # marked terminal nodes idcs
  # Homology classes selection
  react$persistence_selected         <- NULL # selected homology classes idcs
  react$persistence_marked           <- NULL # marked homology classes idcs
  # Trajectories dendrogram
  react$dendrogram                   <- NA   # dendrogram for clustering of trajectories by homology classes
  react$dendrogram_ready             <- FALSE
  react$dendrogram_selection         <- NULL
  react$dendrogram_selected_leaves   <- NULL 
  react$dendrogram_selected_idcs     <- NULL
  react$dendrogram_leaf_perc_cutoff  <- NULL
  react$dendrogram_classes           <- list()
  # 2D trajectories layout
  react$layout_trajectories_flip_colours            <- FALSE
  react$layout_trajectories_highlight_in_background <- FALSE
  # Marker expression and population trackers
  react$trackers_n_segments          <- 20 # number of pseudotime segments in markers tracking
  react$trackers_scaling_exponent    <- 1 # scaling exponent for pseudotime segmentation in markers tracking
  react$tracked_markers.A            <- NULL # names of markers of interest in group A
  react$tracked_markers_stats.A      <- NULL
  react$tracked_markers_stats.B      <- NULL
  react$tracked_markers_parts_to_remove.A <- NULL
  react$tracked_markers_parts_to_remove.B <- NULL
  react$tracked_markers_ready.A      <- FALSE
  react$tracked_markers_ready.B      <- FALSE
  react$tracked_markers_pseudotime_bounds.A <- NULL
  react$tracked_markers_pseudotime_bounds.B <- NULL
  react$tracked_markers_last_removed.A <- NULL
  react$tracked_markers_last_removed.B <- NULL
  react$pseudotime_highlight_bounds  <- NULL
  react$tracked_populations.A        <- NULL
  react$tracked_populations.B        <- NULL
  react$tracked_populations_stats.A  <- NULL
  react$tracked_populations_stats.B  <- NULL
  react$tracked_populations_log2_transform <- FALSE
  # Marked trajectories manipulation & export
  react$trajectories_marked.A        <- NULL
  react$trajectories_marked_leaves.A <- NULL
  react$trajectories_marked.B        <- NULL
  react$trajectories_marked_leaves.B <- NULL
  react$trajectories_junk.A          <- integer(0)
  react$trajectories_junk.B          <- integer(0)
  react$trajectories_group           <- 'A'
  react$output_ff                    <- NULL
  react$trajectories_pinned_batches_count <- 0
  react$trajectories_random_walks    <- NULL
  react$trajectories_to_pin          <- NULL
  react$trajectories_pinned          <- NULL
  
  ## Scale 2D projection
  layout[, 1]         <- layout[, 1] - min(layout[, 1]); layout[, 1] <- layout[, 1] / max(layout[, 1])
  layout[, 2]         <- layout[, 2] - min(layout[, 2]); layout[, 2] <- layout[, 2] / max(layout[, 2])
  layout.df           <- data.frame(layout) # for compatibility with brushOpts
  colnames(layout.df) <- c('X', 'Y')
  
  ## Set up colour palette for plotting annotated populations
  gating_palette       <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
  gating_palette       <- unlist(mapply(RColorBrewer::brewer.pal, gating_palette$maxcolors, rownames(gating_palette)))[-1]
  gating_palette       <- gating_palette[1:length(unique(tv$labels))]
  gating_colour_vector <- gating_palette[as.numeric(as.factor(tv$labels))]
  
  ## Find all trajectories' terminal nodes
  termini              <- c(tv$walks$v[tv$walks$starts[-1] - 1], tail(tv$walks$v, 1))
  termini.unique       <- unique(termini)
  
  ### OUTPUTS & OBSERVERS
  
  ## Help pop-up
  observeEvent(input$btn_help, {
    showModal(modalDialog(
      title = HTML('<h3>How to use the <i>tviblindi</i> Shiny interface</h3>'),
      size = 'l',
      HTML("This interface is a tool to allow for discrimination within a set of canonical developmental trajectories in input data, as well as between canonical and aberrant trajectories. Some trajectories will be biologically relevant, whereas others might not reflect the underlying biology.
             <br><br>
             First, select terminal nodes of simulated random walks in the left pane (tab <i>Terminal nodes selection</i>). This is done by clicking and drawing a rectangular selection. To mark selected points, press the <b>PLUS</b> button underneath. To remove marked points (and start selecting over again), press the <b>FIRE</b> button underneath. Once you are satisfied with your selection (the marked terminal nodes will be clustered together!), press the <b>THUMBS-UP</b> button to compute the relevant triangulation and create trajectory representations. Additionally, press the <b>ENLARGE</b> button in the right part of the left panel to view the dimension reduction layout coloured according to annotated cell populations.
             <br><br>
             Second, select the tab <i>Homology classes by persistence selection</i> in the left pane. Here, select relevant homology classes to be used for hierarchical classification of trajectories. (For our analysis, these correspond roughly to more or less significant holes in the high-dimensional data.) Persistence increases upward. Again, you can make multiple piecewise selections by consecutively drawing rectangles and pressing <b>PLUS</b> to mark the selected points, or you can discard the marked points by pressing the <b>FIRE</b> button.
             <br><br>
             Third, select and mark trajectories of interest using the dendrogram which appears in the middle pane. Each bifurcation in the dendrogram corresponds to difference in path classification with regard to a single homology class. The horizontal coordinate of each bifurcation corresponds to filtration value of addition of a simplex associated with the death of that homology class during filtration. In other words, as we move toward the right, trajectories as clustered together based on the differences in how they circumnavigate around increasingly prominent holes in the high-dimensional space. During the marking of trajectories using the <b>PLUS</b> button, be sure to select the desired category (<b>A</b> or <b>B</b>) for distinguishing between two collections of marked trajectories, as desired.
             <br><br>
             Fourth, inspect the projection of marked trajectories in the 2D layout in the right pane. Category <b>A</b> trajectories are drawn in blue, whereas category <b>B</b> trajectories are drawn in blue. By default, red trajectories are drawn on the top. To flip this ordering, press the <b>FLIP</b> button (with the black-white circular icon) beneath the 2D layout plot.
             <br><br>
             Fifth, inspect the progression of marker expression in either category of trajectories by selecting the <i>Marker expression tracking</i> tab and entering marker(s) of interest in the text boxes below the 2D projection (category <b>A</b> is above category <b>B</b> here). By default, we separate the progression into 20 equally large segments by exponentially scaled pseudotime values. Both the scaling exponent and the number of segments can be adjusted (however, the default settings should give sensible results in most cases). If a single marker of interest is entered, all the trajectories are included in the diagram. In addition, you can remove (un-mark) some trajectories by drawing a selection rectangle in the progression diagram and pressing the <b>LIGHTNING BOLT</b> button. All trajectories passing through the drawn area are 'zapped' (the number of 'zapped' trajectories is printed in the bottom part of the center pane). You cannot zap all pathways in a category. You can undo a single last zap by pressing the <b>BACKWARD</b> button.
             If multiple markers of interest are chosen, the mean value for each segment for each marker is displayed. In that case zapping is not possible. If you are interested in certain pseudotime segments (for example, corresponding to some sudden spike in marker expressions), select points in those segments and press the <b>FLAG</b> button to temporaritly overlay the trajectories layout with a highlight of those points which belong in the corresponding psedotime segments. After pressing the <b>FLAG</b> button, the selection is automatically cleared, therefore clicking the button again will undo the highlighting. Alternatively, press the <b>CROSS</b> button under the 2D trajectories layout to remove the highlighting. In order to switch between a background and foreground positioning of the highlight, press the <b>GLASSES</b> button. Pseudotime segments highlighting works both when viewing single and multiple markers progressions.
             <br><br>
             Sixth, inspect the counts of select annotated cell populations in each pseudotime segment by selecting the <i>Population tracking</i> tab. This is analogical to the marker tracking part. If you wish to display log2 of cell counts, check the <i>log2</i> checkbox.
             <br><br>
             Seventh, export an enhanced FCS file. In the simplest use-case, you can press the <b>SAVE</b> button in the middle pane straightaway (without following any of the above steps) to append two artificial channels to your FCS, containing the 2D layout coordinates displayed in the left and right panes. If you want to append marked trajectories in either of the categories, along with pseudotime values, press the <b>PIN</b> button next to the header for either category. Then, press the <b>SAVE</b> button. To un-pin all the batches of vertices pinned so far, press the <b>GARBARGE</b> button next to it."
    )))
  })
  
  ## Terminal nodes picker
  # Layout
  output$plot_termini <- renderPlot({
    # Compute pseudotime colouring
    psc   <- as.numeric(as.factor(react$pseudotime$res))
    psc   <- psc / max(psc)
    psc   <- psc * 10000 + 1
    cols  <- gplots::greenred(10500)[psc]
    
    par(mar = c(1, 1, 1, 1))
    plot(layout, col = scales::alpha(cols, .05), axes = FALSE, xlab = '', ylab = '', pch = 20, cex = .3, xlim = c(0, 1), ylim = c(0, 1))
    # Plot origin and terminal nodes
    points(layout[tv$origin, 1], layout[tv$origin, 2], col = scales::alpha('yellow', .75), cex = 3, pch = 15)
    if (!is.null(tv$ShowAllFates) && tv$ShowAllFates) {
      if (length(unique(tv$fates)) == 1) {
        points(layout[unique(tv$states), 1], layout[unique(tv$fates), 2], col = scales::alpha('grey', .75), cex = 1.5, pch = 20)
      } else {
        points(layout[unique(tv$fates), ], col = scales::alpha('grey', .75), cex = 1.5, pch = 20)
      }
    }
    if (length(termini.unique) == 1) {
      points(layout[termini.unique, 1], layout[termini.unique, 2], col = scales::alpha('purple', .75), cex = 3, pch = 20) 
    } else {
      points(layout[termini.unique, ], col = scales::alpha('purple', .75), cex = 3, pch = 20)
    }
  })
  
  # Buttons
  observeEvent(input$btn_termini_clear_termini, {
    react$termini_marked          <- data.frame()
    output$log_persistence_marked <- renderText('')
  })
  
  observeEvent(input$btn_termini_update_walks_by_termini, {
    if (length(react$termini_marked) > 0) {
      updated <- .update_walks_by_termini(tv               = tv,
                                          pseudotime       = react$pseudotime,
                                          marked_termini   = react$termini_marked,
                                          termini_per_path = termini)
      react$trajectories_random_walks <- updated$random_walks
      react$representations           <- updated$repre
      react$persistence               <- updated$pers
      react$persistence_diagram       <- updated$pd
      react$dendrogram_ready          <- TRUE
      react$persistence_selected      <- NULL
      react$persistence_marked        <- NULL
      output$log_persistence_marked   <- renderPrint({ cat('Termini were reset\n') })
      react$trajectories_marked.A     <- NULL
      react$trajectories_marked.B     <- NULL
      react$tracked_markers_last_removed.A <- NULL
      react$tracked_markers_last_removed.B <- NULL
    }
  })
  
  # Logs
  output$log_termini_selected <- renderPrint({
    react$termini_selected <- as.numeric(rownames(brushedPoints(layout.df[termini.unique, ], input$selector_termini, xvar = 'X', yvar = 'Y')))
    cat(sapply(react$termini_selected, function(s) paste0(s, ' (', sum(termini == s), ' pathways)')), sep = '\n')
  })
  observeEvent(input$btn_termini_mark_termini, {
    react$termini_marked <- unique(unlist(c(react$termini_marked, react$termini_selected)))
    output$log_termini_marked <- renderPrint({ cat(react$termini_marked, sep = '\n') })
  })
  
  ## Gating layout pop-up
  output$plot_gating_layout <- renderPlot({
    par(mar = c(1, 1, 1, 1))
    gating_layout.colours <- gating_colour_vector
    labels.aligned        <- levels(tv$labels)
    colours.aligned       <- gating_palette[1:length(labels.aligned)]
    # Re-order labels by pseudotime
    average_pseudotime <- c()
    for (label in labels.aligned) {
      idcs        <- which(tv$labels == label)
      idcs.sample <- sample(idcs, min(1000, length(idcs)))
      average_pseudotime <- c(average_pseudotime, mean(tv$pseudotime$res[idcs.sample]))
    }
    ordering <- order(average_pseudotime)
    labels.aligned  <- labels.aligned[ordering]
    colours.aligned <- colours.aligned[ordering]
    
    pch <- rep(20, length(tv$labels)) # symbols (can be set to different for each population here)
    par(xpd = TRUE, mar = c(2, 2, 2, 50))
    
    # Plot ungated events distinctly
    which.ungated <- which(labels.aligned %in% c('ungated', '*ungated*', 'nic'))
    if (length(which.ungated) == 1) {
      label.ungated                  <- labels.aligned[which.ungated]
      idcs.ungated                   <- which(tv$labels == label.ungated)
      colours.aligned[which.ungated] <- gating_colour_vector[idcs.ungated] <- 'black'
      plot(layout[idcs.ungated, ], col = scales::alpha('black', .4), axes = FALSE, xlab = '', ylab = '', pch = '.', cex = .5, xlim = c(0, 1), ylim = c(0, 1))
      points(layout[-idcs.ungated, 1], layout[-idcs.ungated, 2], col = scales::alpha(gating_colour_vector[-idcs.ungated], 0.4), pch = 20, cex = .18)
    } else {
      plot(layout, col = scales::alpha(gating_colour_vector, .35), axes = FALSE, xlab = '', ylab = '', pch = 20, cex = .18, xlim = c(0, 1), ylim = c(0, 1))
    }
    # Plot origin and terminal nodes
    points(layout[tv$origin, 1], layout[tv$origin, 2], col = scales::alpha('yellow', .75), cex = 3, pch = 15)
    if (length(termini.unique) == 1) {
      points(layout[termini.unique, 1], layout[termini.unique, 2], col = scales::alpha('purple', .75), cex = 3, pch = 20) 
    } else {
      points(layout[termini.unique, ], col = scales::alpha('purple', .75), cex = 3, pch = 20)
    }
    # Plot legend for gates colouring
    legend(x = 1.08, y = 1, legend = labels.aligned, fill = colours.aligned, cex = 1.5, bty = 'n')
  })
  
  ## Homology class picker (persistence)
  # Layout
  output$plot_persistence <- renderPlot({
    if (!is.null(react$persistence)) {
      ggplot(react$persistence_diagram, aes(x = xplot, y = yplot, colour = -yplot)) +
        scale_colour_gradientn(colours = rainbow(5)) +
        geom_point(size = .8 + react$persistence_diagram$yplot * 5, alpha = .7) +
        theme_light() + theme(legend.position = 'none') +
        xlab('(Birth + Death) / 2') + ylab('(Death - Birth) / 2')
    }
  })
  
  # Buttons & logs
  output$log_persistence_selected <- renderPrint({
    if (!is.null(react$persistence)) {
      react$persistence_selected <- brushedPoints(react$persistence_diagram, input$selector_persistence)
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
    react$persistence_marked      <- data.frame()
    output$log_persistence_marked <- renderText('No homology classes marked')
    react$trajectories_marked.A   <- NULL
    react$trajectories_marked.B   <- NULL
    session$resetBrush('selector_dendrogram')
  })
  
  ## Trajectories dendrogram
  # Plot
  output$plot_dendrogram <- renderPlot({
    par(mar = c(0, 0, 0, 0))
    if (react$dendrogram_ready) {
      idcs_selected <- as.numeric(unlist(rownames(react$persistence_marked))) # selected homology classes idcs
      perc          <- react$dendrogram_leaf_perc_cutoff
      isolate({
        if (length(idcs_selected) > 0 && !is.null(react$representations)) {
          simplices_selected            <- react$persistence$inds$death[idcs_selected]
          react$representations.reduced <- lapply(react$representations, function(x) x[x %in% simplices_selected])
          dendrogram_plot               <- .trajectories_dendrogram(pers          = react$persistence,
                                                                    repre.reduced = react$representations.reduced,
                                                                    perc          = perc,
                                                                    out.hclust    = react$dendrogram,
                                                                    out.classif   = react$dendrogram_classes)
          
          plot(dendrogram_plot)
          return()
        }
      })
    }
    .draw_placeholder()
  })
  
  # Leaf magnitude cutoff slider (in percentages)
  observeEvent(input$slider_dendrogram_leaf_cutoff, {
    react$dendrogram_leaf_perc_cutoff <- input$slider_dendrogram_leaf_cutoff
  })
  
  # Buttons & logs
  observeEvent(input$btn_dendrogram_mark_leaves, {
    if (react$trajectories_group == 'A') {
      react$trajectories_marked.A        <- unique(unlist(c(react$trajectories_marked.A,        react$dendrogram_selection)))
      react$trajectories_marked_leaves.A <- unique(unlist(c(react$trajectories_marked_leaves.A, react$dendrogram_selected_leaves.A)))
      react$trajcetories_marked_idcs.A   <- unique(unlist(c(react$trajectories_marked_idcs.A),  react$dendrogram_selected_idcs.A))
    } else if (react$trajectories_group == 'B') {
      react$trajectories_marked.B        <- unique(unlist(c(react$trajectories_marked.B,        react$dendrogram_selection)))
      react$trajectories_marked_leaves.B <- unique(unlist(c(react$trajectories_marked_leaves.B, react$dendrogram_selected_leaves.B)))
      react$trajcetories_marked_idcs.B   <- unique(unlist(c(react$trajectories_marked_idcs.B),  react$dendrogram_selected_idcs.B))
    }
  })
  
  observeEvent(input$btn_dendrogram_clear_marked_leaves, {
    if (react$trajectories_group == 'A') {
      react$trajectories_marked.A          <- react$trajectories_marked_leaves.A <- react$trajcetories_marked_idcs.A <- NULL
      react$trajectories_junk.A            <- integer(0)
      react$tracked_markers_last_removed.A <- NULL
    } else if (react$trajectories_group == 'B') {
      react$trajectories_marked.B          <- react$trajectories_marked_leaves.B <- react$trajcetories_marked_idcs.B <- NULL
      react$trajectories_junk.B            <- integer(0)
      react$tracked_markers_last_removed.B <- NULL
    }
  })
  
  observeEvent(react$trajectories_to_pin, {
    if (is.null(react$output_ff)) {
      if (is.null(event_sel)) {
        layout_X <- layout.df$X * 100
        layout_Y <- layout.df$Y * 100
      } else {
        layout_X <- layout_Y <- rep(-100, nrow(input_ff))
        layout_X[event_sel] <- layout.df$X * 100
        layout_Y[event_sel] <- layout.df$Y * 100
      }
      react$output_ff <- fcs.add_col(
        fcs.add_col(
          input_ff, layout_X, colname = 'dimension_reduction_1'
        ),
        layout_Y, colname = 'dimension_reduction_2'
      )
    }
    react$trajectories_pinned_batches_count <- react$trajectories_pinned_batches_count + 1
    walks        <- list()
    walks$v      <- unlist(react$trajectories_random_walks)
    walks$starts <- c(1, cumsum(sapply(react$trajectories_random_walks, length)) + 1)
    react$output_ff <- .add_path_info_to_fcs(ff                       = react$output_ff,
                                             pseudotime               = react$pseudotime,
                                             all_walks                = walks,
                                             trajectories_of_interest = react$trajectories_to_pin,
                                             id                       = react$trajectories_pinned_batches_count,
                                             event_sel                = event_sel)
  })
  
  observeEvent(react$trajectories_pinned_batches_count, {
    count <- react$trajectories_pinned_batches_count
    output$log_pinned_batches_count <- renderPrint({
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
  
  observeEvent(input$btn_trajectories_export_fcs, {
    if (is.null(react$output_ff)) {
      if (is.null(event_sel)) {
        layout_X <- layout.df$X * 100
        layout_Y <- layout.df$Y * 100
      } else {
        layout_X <- layout_Y <- rep(-100, nrow(input_ff))
        layout_X[event_sel] <- layout.df$X * 100
        layout_Y[event_sel] <- layout.df$Y * 100
      }
      react$output_ff <- fcs.add_col(
        fcs.add_col(
          input_ff, layout_X, colname = 'dimension_reduction_1'
        ),
        layout_Y, colname = 'dimension_reduction_2'
      )
    }
    showModal(export_fcs_modal())
  })
  
  observeEvent(input$btn_trajectories_clear_pinned_trajectories, {
    if (is.null(react$output_ff)) {
      if (is.null(event_sel)) {
        layout_X <- layout.df$X * 100
        layout_Y <- layout.df$Y * 100
      } else {
        layout_X <- layout_Y <- rep(-100, nrow(input_ff))
        layout_X[event_sel] <- layout.df$X * 100
        layout_Y[event_sel] <- layout.df$Y * 100
      }
      react$output_ff <- fcs.add_col(
        fcs.add_col(
          input_ff, layout_X, colname = 'dimension_reduction_1'
        ),
        layout_Y, colname = 'dimension_reduction_2'
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
    if (any(!is.na(react$dendrogram))) {
      df  <- data.frame(y = seq(0, 1, length.out = length(react$dendrogram_classes)))
      pts <- brushedPoints(df, input$selector_dendrogram, yvar = 'y')
      pts <- rev(as.numeric(rownames(pts)))

      react$dendrogram_selected_leaves <- names(react$dendrogram_classes)[pts]
      react$dendrogram_selected_idcs   <- as.vector(unlist(react$dendrogram_classes[pts]))
      react$dendrogram_selection       <- as.numeric(unlist(react$dendrogram_classes[pts]))

      sizes <- sapply(react$dendrogram_classes[pts], length)
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
    par(mar = c(0, 0, 0, 0))
    .plot_trajectories(layout,
                       react$trajectories_random_walks,
                       react$trajectories_marked.A,
                       react$trajectories_marked.B,
                       flip_colours = react$layout_trajectories_flip_colours,
                       pseudotime_highlight_bounds = react$pseudotime_highlight_bounds,
                       pseudotime = if (is.null(react$pseudotime_highlight_bounds)) { NULL } else { react$pseudotime },
                       highlight_in_background = react$layout_trajectories_highlight_in_background)
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
    if (input$input_trackers_scaling_exponent > 3 && input$input_trackers_scaling_exponent < 1000) {
      react$trackers_scaling_exponent <- as.integer(input$input_trackers_scaling_exponent)
    }
  })
  observeEvent(input$input_trackers_n_segments, {
    if (input$input_trackers_n_segments > 3 && input$input_trackers_n_segments < 1000) {
      react$trackers_n_segments <- as.integer(input$input_trackers_n_segments)
    }
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
    pts <- brushedPoints(react$tracked_markers_stats.A,
                         input$selector_tracked_markers.A,
                         xvar = 'segment', yvar = 'expression')
    if (nrow(pts) > 0) {
      highlighted_segments <- sort(unique(pts$segment))
      react$pseudotime_highlight_bounds <- react$tracked_markers_pseudotime_bounds.A[c(min(highlighted_segments) - 1, max(highlighted_segments))]
    } else {
      react$pseudotime_highlight_bounds <- NULL
    }
    session$resetBrush('selector_tracked_markers.A')
  })
  
  observeEvent(input$btn_tracked_markers_highlight_segments.B, {
    pts <- brushedPoints(react$tracked_markers_stats.B,
                         input$selector_tracked_markers.B,
                         xvar = 'segment', yvar = 'expression')
    if (nrow(pts) > 0) {
      highlighted_segments <- sort(unique(pts$segment))
      react$pseudotime_highlight_bounds <- react$tracked_markers_pseudotime_bounds.B[c(min(highlighted_segments), max(highlighted_segments))]
    } else {
      react$pseudotime_highlight_bounds <- NULL
    }
    session$resetBrush('selector_tracked_markers.B')
  })
  
  # Plots
  output$plot_tracked_markers.A <- renderPlot({
    if (!is.null(react$trajectories_marked.A) && !is.null(react$tracked_markers.A)) {
      react$tracked_markers_ready.A <- TRUE
      p <- .plot_tracked_markers(react$trajectories_random_walks,
                                 react$trajectories_marked.A,
                                 tv,
                                 react$pseudotime,
                                 markers  = react$tracked_markers.A,
                                 n.part   = react$trackers_n_segments,
                                 exp.part = react$trackers_scaling_exponent)
      react$tracked_markers_stats.A <- p$stats
      react$tracked_markers_pseudotime_bounds.A <- p$pseudotime_bounds
      p$plot
    } else {
      react$tracked_markers_ready.A <- FALSE
    }
  })
  output$plot_tracked_markers.B <- renderPlot({
    if (!is.null(react$trajectories_marked.B) && !is.null(react$tracked_markers.B)) {
      react$tracked_markers_ready.A <- TRUE
      p <- .plot_tracked_markers(react$trajectories_random_walks,
                                 react$trajectories_marked.B,
                                 tv,
                                 react$pseudotime,
                                 markers  = react$tracked_markers.B,
                                 n.part   = react$trackers_n_segments,
                                 exp.part = react$trackers_scaling_exponent)
      react$tracked_markers_stats.B <- p$stats
      react$tracked_markers_pseudotime_bounds.B <- p$pseudotime_bounds
      p$plot
    } else {
      react$tracked_markers_ready.A <- FALSE
    }
  })
  
  ## Population composition trackers
  # Inputs
  observeEvent(input$check_tracked_populations_log2_transform, {
    react$tracked_populations_log2_transform <- input$check_tracked_populations_log2_transform
  })
  observe({
    updateSelectInput(session, 'input_tracked_populations.A', choices = labels.unique)
    updateSelectInput(session, 'input_tracked_populations.B', choices = labels.unique)
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
    if (!is.null(react$trajectories_marked.A) && !is.null(react$tracked_populations.A)) {
      p <- .plot_tracked_populations(react$trajectories_random_walks,
                                     react$trajectories_marked.A,
                                     tv,
                                     react$pseudotime,
                                     populations    = react$tracked_populations.A,
                                     n.part         = react$trackers_n_segments,
                                     exp.part       = react$trackers_scaling_exponent,
                                     log2_transform = react$tracked_populations_log2_transform)
      react$tracked_populations_stats.A <- p$stats
      p$plot
    }
  })
  output$plot_tracked_populations.B <- renderPlot({
    if (!is.null(react$trajectories_marked.B) && !is.null(react$tracked_populations.B)) {
      p <- .plot_tracked_populations(react$trajectories_random_walks,
                                     react$trajectories_marked.B,
                                     tv,
                                     react$pseudotime,
                                     populations    = react$tracked_populations.B,
                                     n.part         = react$trackers_n_segments,
                                     exp.part       = react$trackers_scaling_exponent,
                                     log2_transform = react$tracked_populations_log2_transform)
      react$tracked_populations_stats.B <- p$stats
      p$plot
    }
  })
}
