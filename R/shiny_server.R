shiny_server <- function(  input,
                           output,
                           session  ) {

    shiny_inputs_dir <- "tviblindi_tmp"

    cancel.onSessionEnded <- session$onSessionEnded(function() {
        message(paste0('tviblindi Shiny session ended'))
    })

    message(paste0('Running tviblindi Shiny UI, working directory is ', getwd()))

    tv_name <- readRDS(file.path(shiny_inputs_dir, 'tv.RDS'))
    tv <- get(tv_name, parent.env(environment()))

    layout           <- tv$layout

    event_sel        <- readRDS(file.path(shiny_inputs_dir, 'event_sel.RDS'))

    R                        <- reactiveValues()
    R$pseudotime             <- tv$pseudotime
    R$pers                   <- NULL
    R$repre                  <- NULL
    R$repre.reduced          <- NULL # only include some homology classes
    R$term.selection         <- NULL
    R$term.marked            <- NULL
    R$term.show_gating       <- FALSE
    R$pers.selection         <- NULL
    R$pers.marked            <- data.frame()
    R$markers.n_segments.A   <- NULL
    R$markers.scale_exp.A    <- NULL
    R$markers.n_segments.B   <- NULL
    R$markers.scale_exp.B    <- NULL
    input_ff                 <- flowCore::read.FCS(readRDS(file.path(shiny_inputs_dir, 'input_fcs_path.RDS')))
    R$dendro                 <- NA
    R$dendro_ready           <- FALSE
    R$dendro_data            <- NULL
    R$dendro.selection_nodes <- NULL
    R$dendro.selection_idcs  <- NULL
    R$marked.A               <- NULL
    R$marked_nodes.A         <- NULL
    R$marked_idcs.A          <- NULL
    R$marked.B               <- NULL
    R$marked_nodes.B         <- NULL
    R$marked_idcs.B          <- NULL
    R$junk.A                 <- integer(0)
    R$junk.B                 <- integer(0)
    R$dendro.perc            <- NULL
    R$dendro.classes         <- list()
    R$pathways_categ         <- "A"
    R$lazy_plotting          <- TRUE

    R$output_ff              <- NULL
    R$save_count             <- 0

    R$random_walks           <- NULL

    R$to_append              <- NULL
    R$to_append.A            <- NULL
    R$to_append.B            <- NULL
    markers                  <- colnames(tv$data)
    R$flip_colours           <- FALSE
    R$markers.selected.A     <- NULL
    R$markers.selected.B     <- NULL
    R$expression.stats.A     <- NULL
    R$expression.stats.B     <- NULL
    R$markers.brushed.A      <- NULL
    R$markers.brushed.B      <- NULL

    layout[, 1]      <- layout[, 1] - min(layout[, 1]); layout[, 1] <- layout[, 1] / max(layout[, 1])
    layout[, 2]      <- layout[, 2] - min(layout[, 2]); layout[, 2] <- layout[, 2] / max(layout[, 2])
    layout.df        <- data.frame(layout)
    colnames(layout.df) <- c("X", "Y")

    gating_palette <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
    gating_palette <- unlist(mapply(RColorBrewer::brewer.pal, gating_palette$maxcolors, rownames(gating_palette)))[-1]
    gating_palette <- gating_palette[1:length(unique(tv$labels))]
    gating_colour_vector <- gating_palette[as.numeric(as.factor(tv$labels))]

    term <- c(tv$walks$v[tv$walks$starts[-1] - 1], tail(tv$walks$v, 1))

    ## Interactive persistence diagram
    output$pers_plot <- renderPlot({
        if (!is.null(R$pers)) {
            ggplot(R$pd, aes(x = xplot, y = yplot, color = -yplot)) +
                scale_color_gradientn(colours = rainbow(5)) +
                geom_point(size = 0.8 + R$pd$yplot * 5, alpha = 0.7) +
                theme_light() + theme(legend.position = "none") +
                xlab("(Birth + Death) / 2") + ylab("(Death - Birth) / 2")
        }
    })

    output$term_plot <- renderPlot({
        if (!R$term.show_gating) { # show pseudotime
            psc  <- as.numeric(as.factor(R$pseudotime$res))
            psc  <- psc / max(psc)
            psc  <- psc * 10000 + 1
            col  <- greenred(10500)
            colour_vector <- col[psc]
        } else {
            colour_vector <- gating_colour_vector
        }

        ends <- unique(term)

        #background     <- scattermore(layout, xlim = c(0, 1), ylim = c(0, 1), rgba = col2rgb(col[psc], alpha = TRUE))
        #terminal_nodes <- scattermore(layout[ends, ], xlim = c(0, 1), ylim = c(0, 1), rgba = c(127, 0, 255, 180), cex = 8)
        #plot(background)
        #plot(terminal_nodes, add = TRUE)
        #scattermoreplot(x = layout[, 1], y = layout[, 2], col = col[psc], xlab = '', ylab = '', cex = .5, axes = FALSE)
        #scattermoreplot(x = layout[ends, 1], y = layout[ends, 2], cex = 10, add = TRUE)
        par(mar = c(1, 1, 1, 1))
        plot(layout, col = alpha(colour_vector, 0.05), axes = FALSE, xlab = "", ylab = "", pch = 20, cex = .3, xlim = c(0, 1), ylim = c(0, 1))
        #points(layout[origin, ], col = alpha("purple", 0.75), cex = 1,  pch = 8)
        points(layout[unique(term), ], col = alpha("purple", 0.75), cex = 3, pch = 20)
        points(layout[tv$origin, 1], layout[tv$origin, 2], col = alpha("yellow", 0.75), cex = 3, pch = 15)
    })

    output$term_brush_info <- renderPrint({
        # brushed_area <- input$term_brush
        # if (!is.null(brushed_area)) {
        #     print(brushed_area$xmin)
        #     print(brushed_area$ymin)
        #     print(brushed_area$xmax)
        #     print(brushed_area$ymax)
        # }

        R$term.selection <- as.numeric(rownames(brushedPoints(layout.df[unique(term), ], input$term_brush, xvar = "X", yvar = "Y")))
        out <- sapply(R$term.selection, function(s) paste0(s, " (", sum(term == s), " pathways)"))
        cat(out, sep = "\n")
    })

    ## Button: switch between pseudotime and gating colouring
    # observeEvent(input$term_btn_colouring, {
    #     R$term.show_gating <- !R$term.show_gating
    # })

    # observeEvent(input$term_btn_enlarge, {
    #     showModal(large_plot_modal())
    # })

    ## Button: mark selected persistence diagram points
    observeEvent(input$term_btn_add, {
        R$term.marked <- unique(unlist(c(R$term.marked, R$term.selection)))
        output$term_marked_info <- renderPrint({ cat(R$term.marked, sep = "\n") })
    })

    output$large_plot <- renderPlot({

        ## Messy code, but works for now

        par(mar = c(1, 1, 1, 1))

        large_plot.colours <- gating_colour_vector
        colours.aligned    <- unique(large_plot.colours)
        labels.aligned     <- unique(tv$labels)

        ## Get approximate label ordering
        avg_pseudotime <- c()
        for (l in labels.aligned) {
            idcs           <- which(tv$labels == l)
            idcs.sample    <- sample(idcs, min(100, length(idcs)))
            avg_pseudotime <- c(avg_pseudotime, mean(tv$pseudotime$res[idcs.sample]))
        }
        o <- order(avg_pseudotime)
        labels.aligned  <- labels.aligned[o]
        colours.aligned <- colours.aligned[o]

        pch <- rep(20, length(tv$labels))

        par(xpd = TRUE)
        par(mar = c(2, 2, 2, 50), xpd = TRUE)

        ## If there is an ungated population, plot it differently
        which.ungated <- which(labels.aligned %in% c('ungated', '*ungated*', 'nic'))

        if (length(which.ungated) == 1) {
            label.ungated                  <- labels.aligned[which.ungated]
            idcs.ungated                   <- which(tv$labels == label.ungated)
            colours.aligned[which.ungated] <- gating_colour_vector[idcs.ungated] <- 'black'
            plot(layout[idcs.ungated, ], col = alpha('black', 0.4), axes = FALSE, xlab = '', ylab = '', pch = ".", cex = .5, xlim = c(0, 1), ylim = c(0, 1))
            points(layout[-idcs.ungated, 1], layout[-idcs.ungated, 2], col = alpha(gating_colour_vector[-idcs.ungated], 0.4), pch = 20, cex = .18)
        } else {
            plot(layout, col = alpha(gating_colour_vector, 0.35), axes = FALSE, xlab = '', ylab = '', pch = 20, cex = .18, xlim = c(0, 1), ylim = c(0, 1))
        }

        points(layout[unique(term), ], col = alpha("purple", 0.75), cex = 3, pch = 20)
        points(layout[tv$origin, 1], layout[tv$origin, 2], col = alpha("yellow", 0.75), cex = 3, pch = 15)

        legend(x = 1.08, y = 1, legend = labels.aligned, fill = colours.aligned, cex = 1.5, bty = 'n')
    })

    ## Button: clear marked persistence diagram points
    observeEvent(input$term_btn_clear, {
        R$term.marked           <- data.frame()
        output$pers_marked_info <- renderText("")
    })

    observeEvent(input$term_btn_update, {
        if (length(R$term.marked) > 0) {
            updated <- .update_walks(tv               = tv,
                                     pseudotime       = R$pseudotime,
                                     tt               = R$term.marked, # selected terminal nodes
                                     termini_per_path = term)
            R$random_walks <- updated$walks
            R$repre        <- updated$repre
            R$pers         <- updated$pers
            R$pd           <- updated$pd

            R$dendro_ready <- TRUE

            ## Bump
            # tmp <- R$pers.marked
            # R$pers.marked <- NULL
            # R$pers.marked <- tmp
            # R$marked.A <- NULL
            # R$marked.B <- NULL

            R$pers.marked <- NULL
            R$marked.A <- NULL
            R$marked.B <- NULL
            R$pers.selection <- NULL
        }
    })

    ## Selected homology classes log
    output$pers_brush_info <- renderPrint({
        if (!is.null(R$pers)) {
            R$pers.selection <- brushedPoints(R$pd, input$pers_brush)
            cols.graphical   <- which(colnames(R$pers.selection) %in% c("xplot", "yplot"))
            if (nrow(R$pers.selection) < 16) {
                print(R$pers.selection[, -cols.graphical])
            } else {
                print(R$pers.selection[1:15, -cols.graphical])
                cat("...\n")
            }
        }
    })

    ## Button: mark selected persistence diagram points
    observeEvent(input$pers_btn_add, {
        cols.graphical          <- which(colnames(R$pers.selection) %in% c("xplot", "yplot"))
        R$pers.marked           <- .row_unique(rbind(R$pers.marked, R$pers.selection[, -cols.graphical])) # add newly selected
        R$pers.marked           <- R$pers.marked[order(as.numeric(rownames(R$pers.marked))), ]            # order correctly
        output$pers_marked_info <- renderPrint(
            if (nrow(R$pers.marked) < 16) {
                print(R$pers.marked)
            } else {
                print(R$pers.marked[1:15, ])
                cat("...\n")
            }
        )

        R$marked.A       <- NULL
        R$marked.B       <- NULL
        session$resetBrush("dendro_brush")

        dendro_ready <- (!is.null(R$random_walks) && !is.null(R$pers.marked))
    })

    ## Button: clear marked persistence diagram points
    observeEvent(input$pers_btn_clear, {
        R$pers.marked           <- data.frame()
        output$pers_marked_info <- renderText("No points in persistence diagram marked")

        R$marked.A       <- NULL
        R$marked.B       <- NULL
        session$resetBrush("dendro_brush")
    })

    output$walks_status <- renderText("Please select pathways by terminal nodes first...")
    observeEvent(R$random_walks, {
        if (is.null(R$random_walks)) {
            output$walks_status <- renderText("Please select pathways by terminal nodes first...")
        } else {
            output$walks_status <- renderText("")
        }
    })

    ## Interactive pathways dendrogram
    output$dendro_plot <- renderPlot({
        par(mar = c(0, 0, 0, 0))
        if (R$dendro_ready) {
            idcs.selected <- as.numeric(unlist(rownames(R$pers.marked)))
            perc          <- R$dendro.perc
            isolate({
                if (length(idcs.selected) > 0 && !is.null(R$repre)) {
                    simplices.selected <- R$pers$inds$death[idcs.selected]

                    R$repre.reduced    <- lapply(R$repre, function(x) x[x %in% simplices.selected])

                    p                  <- .pathway_dendrogram(pers               = R$pers,
                                                              repre.reduced      = R$repre.reduced,
                                                              perc_cutoff        = perc,
                                                              out.hclust         = R$dendro,
                                                              out.classification = R$dendro.classes)
                    if(!identical(p, FALSE)) {
                        plot(p)
                        return()
                    }
                }
         })
        }
        .plot_placeholder()
    })

    ## Pathway percentage threshold selector
    observeEvent(input$dendro_perc, { R$dendro.perc <- input$dendro_perc })

    ## Selected pathways log
    output$dendro_brush_info <- renderPrint({
        if (any(!is.na(R$dendro))) {
            df                       <- data.frame(y = seq(0, 1, length.out = length(R$dendro.classes)))
            pts                      <- brushedPoints(df, input$dendro_brush, yvar = "y")
            pts                      <- rev(as.numeric(rownames(pts)))
            R$dendro.selection_nodes <- names(R$dendro.classes)[pts]
            R$dendro.selection_idcs  <- as.vector(unlist(R$dendro.classes[pts]))
            R$dendro.selection       <- as.numeric(unlist(R$dendro.classes[pts]))

            sizes <- sapply(R$dendro.classes[pts], length)
            #cat(R$dendro.selection_nodes, sep = "\n")
            cat(sizes, sep = ", ")
        }
    })

    ## Button: mark selected pathways
    observeEvent(input$dendro_btn_add, {
        if (R$pathways_categ == "A") {
            R$marked.A       <- unlist(c(R$marked.A, R$dendro.selection))
            R$marked_nodes.A <- unlist(c(R$marked_nodes.A, R$dendro.selection_nodes))
            R$marked_idcs.A  <- unlist(c(R$marked_idcs.A, R$dendro.selection_idcs))
        } else if (R$pathways_categ == "B") {
            R$marked.B       <- unlist(c(R$marked.B, R$dendro.selection))
            R$marked_nodes.B <- unlist(c(R$marked_nodes.B, R$dendro.selection_nodes))
            R$marked_idcs.B  <- unlist(c(R$marked_idcs.B, R$dendro.selection_idcs))
        }
    })

    ## Button: clear marked pathways
    observeEvent(input$dendro_btn_clear, {
        if (R$pathways_categ == "A") {
            R$marked.A       <- NULL
            R$marked_nodes.A <- NULL
            R$marked_idcs.A  <- NULL
            R$junk.A         <- integer(0)
        } else if (R$pathways_categ == "B") {
            R$marked.B       <- NULL
            R$marked_nodes.B <- NULL
            R$marked_idcs.B  <- NULL
            R$junk.B         <- integer(0)
        }
    })

    output$dendro_append_info <- renderPrint({ cat("(none)") })

    ## Button: pin marked pathways for addition to output .fcs file
    observeEvent(input$dendro_btn_append.A, {
        R$to_append <- sort(unique(R$marked.A))
    })

    observeEvent(input$dendro_btn_append.B, {
        R$to_append <- sort(unique(R$marked.B))
    })

    observeEvent(R$to_append, {
        if (is.null(R$output_ff)) {
            if (is.null(event_sel)) {
                layoutX <- layout.df$X * 100
                layoutY <- layout.df$Y * 100
            } else {
                layoutX <- layoutY <- rep(-1000, nrow(input_ff))
                layoutX[event_sel] <- layout.df$X * 100
                layoutY[event_sel] <- layout.df$Y * 100
            }
            R$output_ff <- fcs.add_col(
                fcs.add_col(
                    input_ff, layoutX, colname = "dimension_reduction_1"
                ),
                layoutY, colname = "dimension_reduction_2"
            )
        }
        R$save_count <- R$save_count + 1
        w <- list()
        w$v      <- unlist(R$random_walks)
        w$starts <- c(1, cumsum(sapply(R$random_walks, length)) + 1)
        R$output_ff <- addPathInfo2Fcs(R$output_ff,
                                        R$pseudotime,
                                        w,
                                        R$to_append,
                                        ID = R$save_count,
                                        event_sel = event_sel)
    })

    observeEvent(R$save_count, {
        output$save_count_info <- renderPrint(cat(paste0(R$save_count, " ", if (R$save_count == 1) { "batch" } else { "batches" }, " pinned")))
    })

    ## Button: open dialog for saving output .fcs file
    observeEvent(input$dendro_btn_save, {
        if (is.null(R$output_ff)) {
            if (is.null(event_sel)) {
                layoutX <- layout.df$X * 100
                layoutY <- layout.df$Y * 100
            } else {
                layoutX <- layoutY <- rep(-1000, nrow(input_ff))
                layoutX[event_sel] <- layout.df$X * 100
                layoutY[event_sel] <- layout.df$Y * 100
            }
            R$output_ff <- fcs.add_col(
                fcs.add_col(
                    input_ff, layoutX, colname = "dimension_reduction_1"
                ),
                layoutY, colname = "dimension_reduction_2"
            )
        }
        if (!is.null(R$pseudotime) && !is.null(R$random_walks)) showModal(save_modal())
    })

    observeEvent(input$dendro_btn_erase_pinned, {
        if (is.null(R$output_ff)) {
            if (is.null(event_sel)) {
                layoutX <- layout.df$X * 100
                layoutY <- layout.df$Y * 100
            } else {
                layoutX <- layoutY <- rep(-1000, nrow(input_ff))
                layoutX[event_sel] <- layout.df$X * 100
                layoutY[event_sel] <- layout.df$Y * 100
            }
            R$output_ff <- fcs.add_col(
                fcs.add_col(
                    input_ff, layoutX, colname = "dimension_reduction_1"
                ),
                layoutY, colname = "dimension_reduction_2"
            )
        }
        R$save_count <- 0
        R$to_append <- NULL
    })

    observeEvent(input$btn_help, {
        showModal(modalDialog(
            title = "Welcome to tviblindi UI",
            "This interface is a tool to allow for discrimination within a set of canonical developmental trajectories in input data, as well as between canonical and aberrant trajectories. Some trajectories will be biologically relevant, whereas others might not reflect the underlying biology.
             First, select terminal nodes of simulated random walks in the left pane (tab Terminal nodes). This is done by clicking and drawing a rectangular selection. To mark selected points, press the PLUS button underneath. To remove marked points (and start over again), press the FIRE button underneath. Once you're satisfied with you selection (the marked terminal nodes will be clustered together), press the THUMBS-UP button to compute the relevant triangulation and create trajectory representations.
             Second, select the tab Persistence in the left pane. Here, select relevant homology classes to be used for hierarchical classification of trajectories. Persistence increases upward. Again, you can make multiple piecewise selections by consecutively drawing rectangles and pressing PLUS to mark the selected points, or you can discard the marked points by pressing FIRE.
             Third, select and mark trajectories of interest using the dendrogram which appears in the middle pane. Each bifurcation in the dendrogram corresponds to difference in path classification with regard to a single homology class. The horizontal coordinate of each bifurcation corresponds to filtration value associated with addition of a simplex associated with the death of that homology class during filtration. During the marking of trajectories using the PLUS button, be sure to select the desired category (A or B) for distinguishing between two collections of marked trajectories, as desired.
             Fourth, inspect the projection of marked trajectories in the 2-D layout in the right pane. Category A trajectories are drawn in blue, whereas category B trajectories are drawn in blue. By default, red trajectories are drawn on the top. To flip this ordering, press the FLIP button (with the black-white circular icon) beneath the 2-D layout.
             Fifth, inspect the progression of marker expression in either category of trajectories by entering marker(s) of interest the text boxes below the 2-D projection (category A is above category B here). By default, we separate the progression into 20 equally large segments by exponentially scaled pseudotime values. Both the scaling and the number of segments can be adjusted (however, the default settings should give sensible results in most cases). If a single marker of interest is entered, all the trajectories are included in the diagram. In addition, you can remove (un-mark) some trajectories by drawing a selection rectangle in the progression diagram and pressing the LIGHTNING BOLT button. All trajectories passing through the drawn area are 'zapped' (the number of 'zapped' trajectories is printed in the bottom part of the center pane). You cannot zap all pathways in a category.
             If multiple markers of interest are chosen, the mean value for each segment for each marker is displayed.
             Sixth, export an enhanced FCS file. In the simplest case, you can press the SAVE button in the middle pane straightaway to append two artificial channels to your FCS, containing the 2-D layout displayed in the left and right pane. If you want to append marked trajectories in either of the categories, along with pseudotime values, press the PIN button next to the header for either category. Then, press the SAVE button. To un-pin any all batches of vertices pinned so far, press the GARBARGE button next to it."
        ))
    })

    ## Dialog: save output .fcs file
    save_modal <- function(failed = FALSE) {
        modalDialog(
            textInput("output_fcs_name", "Output .fcs file name",
                      placeholder = "", value = "output.FCS"),
            if (failed) div(tags$b("Field is empty", style = "color: red;")),

            footer = tagList(
                modalButton("Cancel"),
                actionButton("save_ok", "Save")
            )
        )
    }

    ## Button: okay saving output .fcs file
    observeEvent(input$save_ok, {
        if (input$output_fcs_name != "") {
            flowCore::write.FCS(R$output_ff, input$output_fcs_name)
            removeModal()
        } else {
            showModal(save_modal(failed = TRUE))
        }
    })

    ## Button switch: pathway category (A vs. B)
    observeEvent(input$dendro_btn_categ, {
        R$pathways_categ <- input$dendro_btn_categ
    })

    observeEvent(input$layout_lazy, {
        R$lazy_plotting <- input$layout_lazy
    })

    ## Data layout projection
    output$layout_plot <- renderPlot({
        par(mar = c(0, 0, 0, 0))
        if (!R$lazy_plotting) {
            plot(layout, cex = 0.05, pch = 20, col = alpha("grey", 0.5), axes = FALSE, xlab = "", ylab = "")
            if (!is.null(R$marked.A) || !is.null(R$marked.B))
                .plot_walks(layout, R$random_walks, R$marked.A, R$marked.B, flip_colours = R$flip_colours)
        } else {
            .plot_lazy_walks(layout, R$random_walks, R$marked.A, R$marked.B, flip_colours = R$flip_colours)
        }
    })

    observeEvent(input$layout_btn_flip_colours, {
        R$flip_colours <- !R$flip_colours
    })

    ## Marker selector choices
    observe({
        updateSelectInput(session, "marker_selector.A", choices = markers)
        updateSelectInput(session, "marker_selector.B", choices = markers)
    })

    ## Marker selectors
    observeEvent(input$marker_selector.A, {
        R$markers.selected.A <- input$marker_selector.A
        if (is.null(R$markers.selected.B)) {
            updateSelectInput(session, "marker_selector.B", selected = input$marker_selector.A)
        }
    })
    observeEvent(input$marker_selector.B, {
        R$markers.selected.B <- input$marker_selector.B
        if (is.null(R$markers.selected.A)) {
            updateSelectInput(session, "marker_selector.A", selected = input$marker_selector.B)
        }
    })

    ## Segmentation inputs
    observeEvent(input$n_segments.A, { if (input$n_segments.A > 3 && input$n_segments.A < 1000) R$markers.n_segments.A <- as.integer(input$n_segments.A) })
    observeEvent(input$n_segments.B, { if (input$n_segments.B > 3 && input$n_segments.B < 1000) R$markers.n_segments.B <- as.integer(input$n_segments.B) })

    ## Scaling exponent inputs
    observeEvent(input$scaling_exponent.A, { if (input$scaling_exponent.A >= 0.1 && input$scaling_exponent.A <= 10) R$markers.scale_exp.A <- input$scaling_exponent.A })
    observeEvent(input$scaling_exponent.B, { if (input$scaling_exponent.B >= 0.1 && input$scaling_exponent.B <= 10) R$markers.scale_exp.B <- input$scaling_exponent.B })

    ## Marker expression plots
    output$expression_plot.A <- renderPlot({
        if (!is.null(R$marked.A) && !is.null(R$markers.selected.A)) {
            p                    <- .plot_expression(R$random_walks, R$marked.A, tv, R$pseudotime,
                                                     markers = R$markers.selected.A,
                                                     n.part = R$markers.n_segments.A,
                                                     exp.part = R$markers.scale_exp.A)
            R$expression.stats.A <- p$stats
            p$plot
        }
    })
    output$expression_plot.B <- renderPlot({
        if (!is.null(R$marked.B) && !is.null(R$markers.selected.B)) {
            p                    <- .plot_expression(R$random_walks, R$marked.B, tv, R$pseudotime,
                                                     markers = R$markers.selected.B,
                                                     n.part = R$markers.n_segments.B,
                                                     exp.part = R$markers.scale_exp.B)
            R$expression.stats.B <- p$stats
            p$plot
        }

    })

    observeEvent(input$expression_zap.A, {
        if (length(R$markers.selected.A) == 1) {
            pts             <- brushedPoints(R$expression.stats.A,
                                             input$expression_brush.A,
                                             xvar = "segment", yvar = "expression")

            if (nrow(pts) < length(R$marked.A)) {
                junk            <- R$marked_idcs.A[unique(pts$walk)]

                R$junk.A        <- na.omit(unique(c(R$junk.A, junk)))
                nonjunk         <- !R$marked_idcs.A %in% junk
                R$marked.A      <- na.omit(R$marked.A[nonjunk])
                R$marked_idcs.A <- na.omit(R$marked_idcs.A[nonjunk])
                session$resetBrush("expression_brush.A")
            }
        }
    })

    observeEvent(input$expression_zap.B, {
        if (length(R$markers.selected.B) == 1) {
            pts             <- brushedPoints(R$expression.stats.B,
                                             input$expression_brush.B,
                                             xvar = "segment", yvar = "expression")

            if (nrow(pts) < length(R$marked.B)) {
                junk            <- R$marked_idcs.B[unique(pts$walk)]

                R$junk.B        <- na.omit(unique(c(R$junk.B, junk)))
                nonjunk         <- !R$marked_idcs.B %in% junk
                R$marked.B      <- na.omit(R$marked.B[nonjunk])
                R$marked_idcs.B <- na.omit(R$marked_idcs.B[nonjunk])
                session$resetBrush("expression_brush.B")
            }
        }
    })

    output$marked_A_info <- output$marked_B_info <- renderPrint(cat("(none)"))

    observeEvent(R$marked.A, {
        R$junk.A <- R$junk.A[!R$junk.A %in% R$marked_idcs.A]

        len <- length(R$junk.A)
        z <- if (len > 0) { paste0(" (", len, " zapped)") } else { "" }
        output$marked_A_counts_info <- renderPrint({ cat("N = ", length(R$marked.A), z, sep = "")})

        if (length(R$marked_idcs.A) > 250) {
            output$marked_A_info <- renderPrint({ cat(paste(R$marked_idcs.A[1:250], sep = ", "), paste0("...and ", length(R$marked_idcs.A) - 250, " more"), sep = "\n") })
        } else if (length(R$marked_idcs.A) > 0) {
            output$marked_A_info <- renderPrint({ cat(R$marked_idcs.A, sep = ", ") })
        } else {
            output$marked_A_info <- renderPrint({ cat("(none)") })
        }
    })

    observeEvent(R$marked.B, {
        R$junk.B <- R$junk.B[!R$junk.B %in% R$marked_idcs.B]

        len <- length(R$junk.B)
        z <- if (len > 0) { paste0(" (", len, " zapped)") } else { "" }
        output$marked_B_counts_info <- renderPrint({ cat("N = ", length(R$marked.B), z, sep = "")})

        if (length(R$marked_idcs.B) > 250) {
            output$marked_B_info <- renderPrint({ cat(paste(R$marked_idcs.B[1:250], sep = ", "), paste0("...and ", length(R$marked_idcs.B) - 250, " more"), sep = "\n") })
        } else if (length(R$marked_idcs.B) > 0) {
            output$marked_B_info <- renderPrint({ cat(R$marked_idcs.B, sep = ", ") })
        } else {
            output$marked_B_info <- renderPrint({ cat("(none)") })
        }
    })
}

## Function: remove duplicit rows from data frame
.row_unique <- function(df) {

    hashes <- apply(df, 1, function(x) paste(x, collapse = " "))
    dupes  <- which(duplicated(hashes))
    if (length(dupes) > 0) return(df[-dupes, ])
    df
}

fcs.add_col <- function(ff, new_col, colname = "label") {

    efcs <- ff@exprs

    N <- nrow(efcs)
    len <- length(new_col)

    if (N != len)  stop(paste0("Number of rows of expression matrix is ", N, ", whereas length of new column is ", len, "."))

    params <- ff@parameters
    pd     <- pData(params)
    cols   <- as.vector(pd$name)
    idcs   <- match(cols, pd$name)

    if (any(is.na(idcs))) stop("Invalid column specifier")

    channel_number     <- ncol(ff) + 1
    channel_id         <- paste0("$P", channel_number)
    channel_name       <- colname
    channel_range      <- max(new_col) + 1
    channel_min        <- min(0, min(new_col) - 1)
    plist              <- matrix(c(channel_name, channel_name, channel_range,
                                   channel_min, channel_range - 1))
    rownames(plist)    <- c("name", "desc", "range", "minRange", "maxRange")
    colnames(plist)    <- c(channel_id)
    pd                 <- rbind(pd, t(plist))
    pData(params)      <- pd
    channel_names      <- colnames(efcs)
    efcs.mod           <- cbind(efcs, new_col)
    colnames(efcs.mod) <- c(channel_names, colname)

    ff.mod             <- flowFrame(efcs.mod, params, description = description(ff))

    keyval                                      <- list()
    keyval[[paste0("$P", channel_number, "B")]] <- "32"
    keyval[[paste0("$P", channel_number, "R")]] <- toString(channel_range)
    keyval[[paste0("$P", channel_number, "E")]] <- "0,0"
    keyval[[paste0("$P", channel_number, "N")]] <- channel_name
    keyval[[paste0("$P", channel_number, "S")]] <- channel_name
    keyword(ff.mod)                             <- keyval

    flowCoreP_Rmax <- paste0("flowCore_$P", channel_number, "Rmax")
    flowCoreP_Rmin <- paste0("flowCore_$P", channel_number, "Rmin")

    description(ff.mod)[flowCoreP_Rmax] <- max(20000, description(ff.mod)$`flowCore_$P1Rmax`)
    description(ff.mod)[flowCoreP_Rmin] <- 0

    return(ff.mod)
}

#### Pathway dendrogram

## Function: remove duplicit rows from data frame
.row_unique <- function(df) {

    hashes <- apply(df, 1, function(x) paste(x, collapse = " "))
    dupes <- which(duplicated(hashes))
    if (length(dupes) > 0) {
        return(df[-dupes, ])
    }
    df
}

## Function: list of numeric vectors to list of strings
.label <- function(l, collapse = " ") {

    lapply(l, function(x) {
        out <- paste(x, collapse = collapse)
        if (out == "") return("nil")
        out
    })
}

## Function: find slot in list containing all elements of given vector
.find <- function(what, where) which(sapply(where, function(x) all(what %in% x)))

## Function: hclust reference -> hclust leaf node (recursively)
.ref_to_vals <- function(ref, where) {

    vals <- vector(mode = "integer")
    vec  <- where[ref, ]
    vals <- c(vals, -vec[vec < 0])
    for (i in vec[vec > 0]) vals <- c(vals, .ref_to_vals(i, where))
    vals
}

## Function: create ggplot2 pathway dendrogram & (optionally) output an hclust object
.pathway_dendrogram <- function(pers,
                                repre.reduced,
                                perc_cutoff,
                                out.hclust         = NULL,
                                out.dendro_data    = NULL,
                                out.classification = NULL) {
    if (perc_cutoff == 100) return(FALSE)

    R              <- repre.reduced; rm(repre.reduced)
    Ru             <- c(unique(unlist(R)))

    threshold      <- ceiling(length(R) / 100 * perc_cutoff)
    H <- lapply(Ru, function(u) {
        idcs <- which(pers$inds$death == u)
        if (length(idcs) == 0) return()
        idcs
    })
    H              <- unlist(H)
    H.vals         <- sort(pers$vals$death[H])
    Ru             <- Ru[order(H)]
    Rl             <- .label(R)
    R.unique_idcs  <- sapply(unique(Rl), function(r) which(Rl == r)[1])
    R.unique       <- R[R.unique_idcs]
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
            sonl                        <- paste(son,  n,  sep = "_")
            sonr                        <- paste(son, -n, sep = "_")

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

    .recurse(pres, 1:length(pres$R), length(pres$Ru), "NULL", "ROOT", "NULL", th = threshold)

    leaves  <- pres$leaves
    h       <- unlist(pres$heights)
    m       <- do.call(rbind, pres$merges)[order(h), ]
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

    ## Resolve leaf ordering
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

    o           <- unlist(o)
    tree        <- list()
    tree$merge  <- merge
    tree$order  <- o
    tree$height <- H.vals[h]
    nodes       <- as.vector(t(m))
    labs        <- nodes[nodes %in% leaves]
    tree$labels <- labs
    class(tree) <- "hclust"

    d           <- as.dendrogram(tree)
    data        <- dendro_data(d, type = "rectangle")
    leaves      <- labs

    labs              <- as.character(data$labels$label)
    data$labels$label <- as.character(sapply(as.character(data$labels$label), function(leaf) get(leaf, pres)[[2]]))
    cl                <- lapply(labs, function(leaf) get(leaf, pres)[[3]]) # walk indices per leaf
    names(cl)         <- labs

    if (!is.null(out.classification)) eval.parent(substitute(out.classification <- cl))
    if (!is.null(out.hclust))         eval.parent(substitute(out.hclust         <- tree))
    if (!is.null(out.dendro_data))    eval.parent(substitute(out.dendro_data    <- data))

    p <- ggplot(data$segments) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_text(data = data$labels, aes(x, y, label = label), hjust = 1, angle = 0, size = 3.4) +
        #ylim(min(tree$height) -0.3, max(tree$height)) +
        cowplot::theme_nothing() +
        #labs(x = NULL, y = NULL) + #scale_x_continuous(expand = c(.01, .01)) +
        coord_flip()
}

## Functions: add selected random walks to plot of a projection layout
.plot_walks <- function(X, walks, walk_idcs.A, walk_idcs.B, flip_colours, ...) {

    cols <- c("blue", "darkred")

    j      <- 0
    sel1   <- walks[walk_idcs.A]
    sel2   <- walks[walk_idcs.B]
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
        lines(pts, lty = 1, col = alpha("darkgrey", 0.02), ...)
        points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = alpha("darkgrey", 0.02))
    }

    if (length(sel1) > 0) {
        for (i in sel1) {
            j   <- j + 1
            pts <- X[i, ]
            lines(pts, lty = 1, col = alpha(cols[1], 0.3), ...)
            points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = alpha(cols[1], 0.3))
        }
    }
    if (length(sel1) > 0) {
        for (i in sel2) {
            j   <- j + 1
            pts <- X[i, ]
            lines(pts, lty = 1, col = alpha(cols[2], 0.3), ...)
            points(pts[, 1], pts[, 2], pch = 20, cex = 0.2, col = alpha(cols[2], 0.3))
        }
    }
}

.plot_lazy_walks <- function(X, walks, walk_idcs.A, walk_idcs.B, flip_colours, ...) {

    col1 <- c(34, 87, 201, 255)
    col2 <- c(194, 45, 55, 255)

    # X[, 1] <- X[, 1] - min(X[, 1]); X[, 1] <- X[, 1] / max(X[, 1])
    # X[, 2] <- X[, 2] - min(X[, 2]); X[, 2] <- X[, 2] / max(X[, 2])
    plot(scattermore(X, rgba = c(200, 200, 200, 150)))

    j      <- 0
    sel1   <- walks[walk_idcs.A]
    sel2   <- walks[walk_idcs.B]

    if (flip_colours) {
        tmp <- sel1
        sel1 <- sel2
        sel2 <- tmp
        tmp <- col1
        col1 <- col2
        col2 <- tmp
    }

    if (length(sel1) > 0) {
        pts1 <- vector(mode = "list", length = length(sel1))
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
        pts2 <- vector(mode = "list", length = length(sel2))
        for (idx in 1:length(sel2)) {
            s   <- sel2[[idx]]
            j   <- j + 1
            pts <- X[s, ]

            pts2[[idx]] <- do.call(rbind, interpolate_trajectories(lapply(1:nrow(pts), function(x) pts[x, ])))
        }
        pts2 <- do.call(rbind, pts2)
        plot(scattermore(pts2, rgba = col2, xlim = c(0, 1), ylim = c(0, 1)), add = TRUE)
    }
}

## Functions: plot expression markers across pseudotime for selected walks
.plot_expression <- function(walks,
                             walk_idcs,
                             tv,
                             pseudotime,
                             markers,
                             breaks = NULL,
                             n.part = 10,
                             exp.part = 1,
                             asinh.denominator = 5,
                             show_legend = FALSE) {

    if (length(markers) == 1) {
        return(do.call(.plot_expression.single,
                       as.list(environment())))
    } else {
        return(do.call(.plot_expression.multiple,
                       as.list(environment())))
    }
}

.plot_expression.single <- function(walks,
                                    walk_idcs,
                                    tv,
                                    pseudotime,
                                    markers,
                                    breaks = NULL,
                                    n.part = 10,
                                    exp.part = 1,
                                    asinh.denominator,
                                    show_legend = FALSE) {
    pseudotime <- pseudotime$res
    pseudotime <- pseudotime/max(pseudotime)

    walks      <- walks[walk_idcs]

    progress   <- lapply(walks, function(pts) pseudotime[pts])
    pp         <- unique(unlist(progress))

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
            b <- sapply(0:breaks, function(i) 1 / breaks * i)
            N <- length(b) - 1
        }
    }

    categs <- as.numeric(cut(unlist(progress), breaks = b, include.lowest = TRUE))

    coords <- tv$data[, markers]

    stats  <- lapply(1:N, function(i) {
        inds <- categs == i # pick points on paths by pseudotime increment

        if (!any(inds)) return(NULL)

        which.walks <- unlist(lapply(1:length(walks), function(j) rep(j, length(walks[[j]]))))[inds]
        pts         <- unlist(walks)[unlist(inds)]
        means       <- lapply(sort(unique(which.walks)), function(j) mean(coords[pts[which.walks == j]]))

        v           <- cbind(rep(i, length(means)), sort(unique(which.walks)), unlist(means))
        if (is.null(dim(v))) { names(v) <- c("segment", "walk", "expression") } else { colnames(v) <- c("segment", "walk", "expression") }
        as.data.frame(v)
    })

    stats            <- do.call(rbind, stats)
    stats$walk       <- as.factor(stats$walk)
    stats$segment    <- as.numeric(stats$segment)
    stats$expression <- as.numeric(stats$expression)

    g <- ggplot(data = stats, aes(x = segment, y = expression, group = walk, color = walk)) +
        geom_line(linetype = "dashed") + geom_point() +
        ggtitle(paste0(markers, " expression per walk: means per segment")) +
        labs(subtitle = "Segmented by pseudotime values.") +
        theme_light() +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 20), legend.title = element_text(size = 24), legend.text = element_text(size = 24))
    if (!show_legend) {
        g <- g + theme(legend.position = "none")
    }
    list(plot  = g,
         stats = stats)
}

.plot_expression.multiple <- function(walks,
                                      walk_idcs,
                                      tv,
                                      pseudotime,
                                      markers,
                                      breaks = NULL,
                                      n.part = 10,
                                      exp.part = 1,
                                      asinh.denominator = 5,
                                      show_legend = TRUE) {
    pseudotime <- pseudotime$res
    pseudotime <- pseudotime/max(pseudotime)

    walks    <- walks[walk_idcs]
    progress <- lapply(walks, function(pts) pseudotime[pts])
    pp       <- unique(unlist(progress))

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
            b <- sapply(0:breaks, function(i) 1 / breaks * i)
            N <- length(b) - 1
        }
    }

    categs <- as.numeric(cut(unlist(progress), breaks = b, include.lowest = TRUE))

    coords <- tv$data[, markers]

    stats  <- lapply(markers, function(m) {
        lapply(1:N, function(i) {

            inds <- categs == i # pick points on paths by pseudotime increment
            if (!any(inds)) return(NULL)

            which.walks <- unlist(lapply(1:length(walks), function(j) rep(j, length(walks[[j]]))))[inds]
            pts         <- unlist(walks)[unlist(inds)]
            means       <- sapply(sort(unique(which.walks)), function(j) mean(coords[pts[which.walks == j], m])) %>% mean


            v <- cbind(rep(m, length(means)), rep(i, length(means)), unlist(means))
            if (is.null(dim(v))) { names(v) <- c("marker", "segment", "expression") } else { colnames(v) <- c("marker", "segment", "expression") }
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
        ggtitle(paste0("Multiple markers expression")) +
        labs(subtitle = "Segmented by pseudotime values. ") +
        theme_light() +
        theme(axis.text.x = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 20), legend.title = element_text(size = 24), legend.text = element_text(size = 24))

    list(plot  = g,
         stats = stats)
}

.update_walks <- function(tv, pseudotime, tt, termini_per_path) {

    ## Identify chosen walks
    idcs <- which(termini_per_path %in% tt)

    ## Remove other walks
    walks.selected <- lapply(idcs, function(idx) select_paths_points(tv$walks, idx))
    lens           <- sapply(walks.selected, length)
    walks.selected <- list(v      = unlist(walks.selected),
                           starts = c(1, 1 + cumsum(lens[-length(lens)])))

    ## Put terminal node with highest pseudotime at end of each walk
    if (length(tt) > 1) {
        max_t             <- tt[which.max(pseudotime$res[tt])][1]
        ends              <- c(walks.selected$starts[-1] - 1, length(walks.selected$v))
        walks.selected$v[ends] <- max_t
    }

    ## Cluster and triangulate walks
    withProgress(message = "Triangulating paths", style = "old", expr = {
        walks_clusters <- remove_cycles(contract_walks(walks.selected, tv$clusters), verbose = FALSE)
        triangulation  <- triangulate_pathways(walks_clusters, tv$codes, tv$filtration$cmplx)
    })

    N <- length(idcs)
    repre      <- list()[1:N]
    repre[[1]] <- integer(0)
    j          <- 1

    N <- length(idcs)
    tick <- 1 / N

    withProgress(message = "Computing representations", value = tick, style = "old", expr = {
        for (idx in 2:N) {
            j          <- j + 1
            cycle      <- path_sum(triangulation[[1]], triangulation[[idx]])
            this_repre <- get_rep_straight(cycle, tv$reduced_boundary, tv$boundary, update = TRUE)
            repre[[j]] <- this_repre
            incProgress(tick)
        }
    })

    walks <- lapply(1:N, function(idx) { select_paths_points(walks.selected, idx) })

    pers <- pers_diagram(dBr = tv$reduced_boundary, plot = FALSE, repre = repre)
    pd   <- data.frame(Dimension    = pers$vals$dim,
                       Birth        = pers$vals$birth,
                       Death        = pers$vals$death,
                       BirthSimplex = pers$inds$birth,
                       DeathSimplex = pers$inds$death,
                       xplot        = (pers$vals$birth + pers$vals$death) / 2,
                       yplot        = (pers$vals$death - pers$vals$birth) / 2)
    pd   <- pd[pd$Dimension > 0, ]

    list(walks  = walks,
         repre  = repre,
         pers   = pers,
         pd     = pd)
}

.plot_placeholder <- function() {
    j   <- jpeg::readJPEG(system.file("tree.jpg", package = "tviblindi"), native = TRUE)
    plot(0:1, 0:1, type = "n", ann = FALSE, axes = FALSE)
    rasterImage(j, 0, 0, 1, 1)
}

addPathInfo2Fcs <- function(fcs, pseudotime, walks, walks.selection, event_sel = NULL, ID = NULL) {
    pp            <- unique(select_paths_points(walks, walks.selection))
    out           <- matrix(-1, nrow = nrow(fcs@exprs), ncol = 2)
    colnames(out) <- c(paste(ID, "which_event", sep = "_"), paste(ID, "local_pseudotime", sep = "_"))
    if (is.null(event_sel)) event_sel <- 1:nrow(fcs@exprs)

    out[event_sel[pp], 1] <- 1000
    out[event_sel[pp], 2] <- as.numeric(as.factor(pseudotime$res[pp]))

    make_valid_fcs(cbind(fcs@exprs, out), desc1 = as.character(fcs@parameters@data$desc))
}
