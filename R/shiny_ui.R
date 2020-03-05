#### Shiny interface for tviblindi: UI layout

shiny_ui <- fluidPage(
  ## CSS tags
  # STYLE: tviblindi title font
  tags$style(HTML('#app_title {
                    font-size: 32px;
                    font-family: monospace
                  }')),
  tags$style(type = 'text/css',
             # STYLE: leaf magnitude cutoff slider
             '.irs-slider {
                width: 15px;
                height: 15px;
                top: 22px
             }'),
  tags$head(
    tags$style('#modal_layout_gating .modal-lg {
          width: 1850px
       }
  ')),
  
  tags$head(
    tags$style('
        .shiny-notification {
            position: fixed;
            top: 20%;
            left: 10%;
            font-size: 24px
  ')),
  
  ## MODALS
  shinyBS::bsModal(
    id      = 'modal_layout_gating',
    size    = 'large',
    title   = '2-dimensional layout with annotated populations',
    trigger = 'btn_left_show_gating',
    plotOutput('plot_gating_layout',
               width  = '1800px',
               height = '950px')
  ),
  shinyBS::bsModal(
    id      = 'modal_tracked_markers.A',
    size    = 'large',
    title   = 'Tracking marker expression',
    trigger = 'btn_tracked_markers_enlarge.A',
    plotOutput('plot_tracked_markers_enlarged.A',
               width  = '1000px',
               height = '750px')
  ),
  shinyBS::bsModal(
    id      = 'modal_tracked_markers.B',
    size    = 'large',
    title   = 'Tracking marker expression',
    trigger = 'btn_tracked_markers_enlarge.B',
    plotOutput('plot_tracked_markers_enlarged.B',
               width  = '1000px',
               height = '750px')
  ),
  
  fluidRow(
    ## APP TITLE
    column(
      width = 11,
      titlePanel(
        h3(id = 'app_title', HTML('tvi<b>blindi</b>')),
        windowTitle = 'tviblindi'
      )
    ),
    ## HELP BUTTON
    column(
      width = 1,
      style = 'text-align:  right;
               padding-right: 5px;
               margin-top:   15px',
      actionButton('btn_help', '', icon = icon('question-circle', lib = 'font-awesome'))
    )
  ),
  br(),
  mainPanel(
    style = 'width: 100%',
    ## LEFT PANEL
    column(
      width = 4,
      tabsetPanel(
        ## LEFT PANEL TAB: TERMINAL NODES SELECTOR
        tabPanel(
          'Terminal nodes selection',
          fluidRow(
            plotOutput(
              'plot_termini',
              height = 500,
              brush = brushOpts(
                id     = 'selector_termini',
                fill   = 'yellow',
                stroke = 'yellow'
              )
            )
          ),
          column(
            width = 8,
            style = 'text-aligh:    left;
                     pading-bottom: 10px',
            actionButton('btn_termini_mark_termini',            '', icon = icon('glyphicon glyphicon-plus',      lib = 'glyphicon')),
            actionButton('btn_termini_clear_termini',           '', icon = icon('glyphicon glyphicon-fire',      lib = 'glyphicon')),
            HTML("&nbsp;&nbsp;&nbsp;&nbsp;"),
            actionButton('btn_termini_update_walks_by_termini', '', icon = icon('glyphicon glyphicon-thumbs-up', lib = 'glyphicon'))
          ),
          column(
            width = 4,
            style = 'text-align:    right;
                     padding-bottom: 10px',
            actionButton('btn_left_show_gating',             '', icon = icon('glyphicon glyphicon-fullscreen', lib = 'glyphicon'))
          ),
          fluidRow(
            h4('Selected terminal nodes'),
            verbatimTextOutput('log_termini_selected', placeholder = TRUE),
          ),
          fluidRow(
            h4('Terminal nodes marked for further analysis'),
            verbatimTextOutput('log_termini_marked', placeholder = TRUE)
          )
        ),
        ## LEFT PANEL TAB: PERSISTENCE SELECTOR (HOMOLOGY CLASSES)
        tabPanel(
          'Homology classes by persistence selection',
          fluidRow(
            plotOutput(
              'plot_persistence',
              height = 500,
              brush  = brushOpts(
                id     = 'selector_persistence',
                fill   = 'yellow',
                stroke = 'yellow'
              )
            )
          ),
          actionButton('btn_persistence_mark_classes',  '', icon = icon('glyphicon glyphicon-plus', lib = 'glyphicon')),
          actionButton('btn_persistence_clear_classes', '', icon = icon("glyphicon glyphicon-fire", lib = 'glyphicon')),
          HTML('&nbsp;&nbsp;'),
          textOutput('log_persistence_available', inline = TRUE),
          fluidRow(
            h4('Selected homology classes'),
            verbatimTextOutput('log_persistence_selected', placeholder = TRUE)
          ),
          fluidRow(
            h4('Marked homology classes'),
            verbatimTextOutput('log_persistence_marked', placeholder = TRUE)
          )
        )
      )
    ),
    ## MIDDLE PANEL: TRAJECTORY DENDROGRAM
    column(
      width = 4,
      style = 'padding-left: 50px',
      # Dendrogram of trajectories clustered by marked homology classes
      fluidRow(
        plotOutput('plot_dendrogram',
                   height = 800,
                   brush  = brushOpts(
                     id        = 'selector_dendrogram',
                     fill      = 'yellow',
                     stroke    = 'yellow',
                     direction = 'y'
                   ))
      ),
      # Leaf magnitude cutoff slider (in percentages)
      sliderInput(
        'slider_dendrogram_leaf_cutoff',
        'Min trajectory count % per leaf',
        min = 0, max = 100,
        value = 10.0,
        step  = 1
      ),
      # Trajectory group switch (A vs. B)
      radioGroupButtons(
        'btn_trajectories_group',
        label = NULL,
        choiceValues = c('A', 'B'),
        choiceNames  = c('A', 'B'),
        selected = 'A', status = 'default',
        size = 'sm', direction = 'horizontal',
        justified = TRUE, individual = F
      ),
      tags$script("$(\"input:radio[name='btn_trajectories_group'][value='A']\").parent().css('background-color', '#c2dfff');"),
      tags$script("$(\"input:radio[name='btn_trajectories_group'][value='B']\").parent().css('background-color', '#ffb5c9');"),
      column(
        width = 6,
        actionButton('btn_dendrogram_mark_leaves',         '', icon = icon('glyphicon glyphicon-plus',  lib = 'glyphicon')),
        actionButton('btn_dendrogram_clear_marked_leaves', '', icon = icon('glyphicon glyphicon-fire',  lib = 'glyphicon')),
        HTML('&nbsp;&nbsp;&nbsp;&nbsp;'),
        actionButton('btn_trajectories_export_fcs',          '', icon = icon('glyphicon glyphicon-save',  lib = 'glyphicon')),
        actionButton('btn_trajectories_clear_pinned_trajectories', '', icon = icon('glyphicon glyphicon-trash', lib = 'glyphicon'))
      ),
      column(
        width = 6,
        style = 'text-align:  right;
                 padding-right: 5px;
                 margin-top:   15px',
        textOutput('log_pinned_batches_count')
      ),
      br(), br(),
      fluidRow(
        width = 12,
        h4('Selected dendrogram nodes by counts'),
        verbatimTextOutput('log_dendrogram_selected', placeholder = TRUE)
      ),
      fluidRow(
        style = 'padding-left:  0px;
                 padding-right: 0px',
        column(
          width = 6,
          style = 'padding-left:  0px;
                   padding-right: 5px',
          h4('Marked pathways in A'),
          column(
            width = 6,
            actionButton('btn_trajectories_pin_trajectories.A', '', icon = icon('glyphicon glyphicon-pushpin', lib = 'glyphicon'))
          ),
          column(
            width = 6,
            textOutput('log_trajectories_marked_count.A')
          ),
          br(), br(),
          verbatimTextOutput('log_trajectories_marked.A', placeholder = TRUE)
        ),
        column(
          width = 6,
          style = 'padding-left:  0px;
                   padding-right: 5px',
          h4('Marked pathways in B'),
          column(
            width = 6,
            actionButton('btn_trajectories_pin_trajectories.B', '', icon = icon('glyphicon glyphicon-pushpin', lib = 'glyphicon'))
          ),
          column(
            width = 6,
            textOutput('log_trajectories_marked_count.B')
          ),
          br(), br(),
          verbatimTextOutput('log_trajectories_marked.B', placeholder = TRUE)
        )
      ),
      br(), br(), br(),
    ),
    ## RIGHT PANEL
    column(
      width = 4,
      style = 'padding-left: 50px',
      ## RIGHT PANEL: TRAJECTORIES LAYOUT
      fluidRow(
        plotOutput(
          'plot_layout_trajectories',
          height = 500
        )
      ),
      fluidRow(
        column(
          width = 12,
          style = 'text-align:  right;
                   padding-right: 5px;
                   margin-top:   15px',
          actionButton('btn_layout_trajectories_remove_highlight', label = '', icon = icon('glyphicon glyphicon-remove', lib = 'glyphicon')),
          actionButton('btn_layout_trajectories_highlight_in_background', label = '', icon = icon('glyphicon glyphicon-sunglasses ', lib = 'glyphicon')),
          HTML("&nbsp;&nbsp;"),
          actionButton('btn_layout_trajectories_flip_colours', label = '', icon = icon('glyphicon glyphicon-adjust', lib = 'glyphicon'))
        ),
        ### RIGHT PANEL: TRACKERS
        hr(),
        fluidRow(
          column(
            width = 2,
            style = 'text-align:  right;
                     padding-left: 5px;
                     margin-top:   15px',
            numericInput('input_trackers_scaling_exponent', label = 'Scale exp', min = .1, max = 10, value = 1, step = .1),
          ),
          column(
            width = 2,
            style = 'text-align:  right;
                     padding-left: 5px;
                     margin-top:   15px',
            numericInput('input_trackers_n_segments', label = '# segs', min = 3, max = 1000, value = 20, step = 1)
          ),
          column(
            width = 8,
            style = 'text-align:  right;
                     padding-left: 5px;
                     margin-top:   35px',
            checkboxInput('check_trackers_large_base_size', label = 'Larger text', value = FALSE)
          )
        ),
        tabsetPanel(
          tabPanel(
            title = 'Marker expression tracking',
            br(),
            br(),
            ## RIGHT PANEL TAB: MARKER EXPRESSION TRACKERS
            ## A
            fluidRow(
              column(
                width = 9,
                style = 'padding-left: 5px',
                selectInput('input_tracked_markers.A', label = 'Group A markers of interest', choices = c(), multiple = TRUE, width = '100%')
              ),
              column(
                width = 3,
                style = 'padding-left:2px;
                     margin-top: 26px;
                    text-align: right',
                actionButton('btn_tracked_markers_remove_trajectories.A', label = '',      icon = icon('glyphicon glyphicon-flash',          lib = 'glyphicon')),
                actionButton('btn_tracked_markers_undo_remove_trajectories.A', label = '', icon = icon('glyphicon glyphicon-step-backward',  lib = 'glyphicon')),
                actionButton('btn_tracked_markers_highlight_segments.A',  label = '',      icon = icon('glyphicon glyphicon-flag',           lib = 'glyphicon'))
              )
            ),
            fluidRow(
              plotOutput(
                'plot_tracked_markers.A',
                height = 500,
                brush = brushOpts(
                  id     = 'selector_tracked_markers.A',
                  fill   = 'red',
                  stroke = 'red'
                )
              )
            ),
            ## B
            fluidRow(
              column(
                width = 9,
                style = 'padding-left: 5px',
                selectInput('input_tracked_markers.B', label = 'Group B markers of interest', choices = c(), multiple = TRUE, width = '100%')
              ),
              column(
                width = 3,
                style = 'padding-left: 2px;
                   margin-top:  26px;
                   text-align: right',
                actionButton('btn_tracked_markers_remove_trajectories.B', label = '',      icon = icon('glyphicon glyphicon-flash',          lib = 'glyphicon')),
                actionButton('btn_tracked_markers_undo_remove_trajectories.B', label = '', icon = icon('glyphicon glyphicon-step-backward',  lib = 'glyphicon')),
                actionButton('btn_tracked_markers_highlight_segments.B', label = '',       icon = icon('glyphicon glyphicon-flag',           lib = 'glyphicon'))
              )
            ),
            fluidRow(
              plotOutput(
                'plot_tracked_markers.B',
                height = 500,
                brush = brushOpts(
                  id     = 'selector_tracked_markers.B',
                  fill   = 'red',
                  stroke = 'red'
                )
              )
            )
          ),
          br(),
          tabPanel(
            title = 'Population tracking',
            br(),
            ## RIGHT PANEL TAB: ANNOTATED POPULATIONS COMPOSITION TRACKERS
            ## A
            fluidRow(
              column(
                width = 12,
                checkboxInput('check_tracked_populations_log2_transform', label = 'log2', value = FALSE)
              )
            ),
            fluidRow(
              column(
                width = 12,
                style = 'padding-left: 5px;
                         padding-right: 5px',
                selectInput('input_tracked_populations.A', label = 'Group A populations of interest', choices = c(), multiple = TRUE, width = '100%')
              )
            ),
            fluidRow(
              plotOutput(
                'plot_tracked_populations.A',
                height = 500,
              )
            ),
            ## B
            fluidRow(
              column(
                width = 12,
                style = 'padding-left: 5px;
                         padding-right: 5px',
                selectInput('input_tracked_populations.B', label = 'Group B populations of interest', choices = c(), multiple = TRUE, width = '100%')
              )
            ),
            fluidRow(
              plotOutput(
                'plot_tracked_populations.B',
                height = 500,
              )
            )
          )
        )
    )
  )
))
