## Shiny interface module
## UI layout

require(shiny)
require(shinyWidgets)

shiny_ui <- fluidPage(
  title = 'tviblindi',
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
    tags$script(HTML('
        $(window).resize(function(event){
          var w = $(this).width();
          var h = $(this).height();
          var obj = {width: w, height: h};
          Shiny.onInputChange("windowSizeChange", obj);
        });
      '))
  ),
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
  
  tags$head(
    tags$style(HTML('#btn_left_show_gating {background-color:#337ab7; color:white}')),
    tags$style(HTML('#btn_help {background-color:#337ab7; color:white}')),
  ),
  
  ## MODALS
  shinyBS::bsModal(
    id      = 'modal_layout_gating',
    size    = 'large',
    title   = '2-dimensional layout with annotated populations',
    trigger = 'btn_left_show_gating',
    shinycssloaders::withSpinner(
      plotOutput('plot_gating_layout',
                 width  = '1800px',
                 height = '950px'),
      color='#1d2c8f'
    )
  ),
  fluidRow(
    style = "background:linear-gradient(0deg, rgba(255,255,255,1) 0%, #d6d6d6 100%)",
    ## APP TITLE
    column(
      width = 8,
      column(
        align = 'center',
        style = 'color:white; padding-top:15px; background-color:#337ab7',
        width = 2,
        titlePanel(
          h2(id = 'app_title', HTML('tvi<b>blindi</b>')),
          windowTitle = 'tviblindi'
        )
      ),
      column(
        width = 10,
        style = 'padding-top: 10px; color: #000052',
        fluidRow(
          column(
            width = 3,
            style = 'text-align:    left; padding-left:45px',
            selectInput(
              'input_dimred_method', label = 'LAYOUT TYPE', c('default'),
              selected = NULL,
              multiple = FALSE
            )
          ),
          column(
            width = 3,
            style = 'text-align:    left; padding-left:45px',
            selectInput(
              'input_labels_name', label = 'LABELS VECTOR', c('default'),
              selected = NULL,
              multiple = FALSE
            )
          ),
          column(
            width = 2,
            style = 'text-align:    left; padding-left:40px',
            radioGroupButtons(
              inputId = 'btn_layout_pointsize',
              label = 'POINT SIZE',
              choiceValues = c(.065, .4, 1.1),
              choiceNames = c('S', 'M', 'L'),
              status = 'primary'
            )
          ),
          column(
            width = 1,
            style = 'text-align: left; margin-top: 26px',
            actionButton('btn_left_show_gating', '', icon = icon('glyphicon glyphicon-fullscreen', lib = 'glyphicon'))
          ),
          column(
            width = 3,
            style = 'text-align: left',
            radioGroupButtons(
              'btn_image_export_format',
              label = 'IMAGE EXPORT',
              choiceValues = c('PNG', 'SVG'),
              choiceNames  = c('PNG', 'SVG'),
              selected = 'PNG', status = 'primary',
              size = 's', direction = 'horizontal',
              justified = TRUE, individual = FALSE,
              width = 150
            )
          )
        )
      )  
    ),
    ## ANALYSIS NAME & HELP BUTTON
    span(
      style = 'color:#000052',
      column(
        width = 4,
        style = 'text-align:  right;
               padding-right: 25px;
               margin-top:   26px',
        textOutput(outputId = 'text_analysis_name', inline = TRUE),
        HTML('&nbsp;&nbsp;&nbsp;'),
        actionButton('btn_help', '', icon = icon('question-circle', lib = 'font-awesome'))
      )
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
            shinycssloaders::withSpinner(
              plotOutput(
                'plot_termini',
                height = 500,
                brush = brushOpts(
                  id     = 'selector_termini',
                  fill   = '#0242f2',
                  stroke = '#0242f2'
                )
              ),
              color='#1d2c8f'
            )
          ),
          column(
            width = 4,
            style = 'text-align:    left;
                     pading-bottom: 10px',
            actionButton('btn_termini_mark_termini',            '', icon = icon('glyphicon glyphicon-plus',      lib = 'glyphicon')),
            actionButton('btn_termini_clear_termini',           '', icon = icon('glyphicon glyphicon-fire',      lib = 'glyphicon')),
            HTML("&nbsp;&nbsp;&nbsp;&nbsp;"),
            actionButton('btn_termini_update_walks_by_termini', '', icon = icon('glyphicon glyphicon-thumbs-up', lib = 'glyphicon'))
          ),
          column(
            width = 8,
            style = 'text-align:    right;
                     padding-bottom: 10px',
            actionButton('btn_termini_export_image',           '', icon = icon('glyphicon glyphicon-picture',  lib = 'glyphicon'))
          ),
          fluidRow(
            h4('Selected terminal nodes'),
            verbatimTextOutput('log_termini_selected', placeholder = TRUE)
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
            shinycssloaders::withSpinner(
              plotOutput(
                'plot_persistence',
                height = 500,
                brush  = brushOpts(
                  id     = 'selector_persistence',
                  fill   = 'yellow',
                  stroke = 'yellow'
                )
              ),
              color='#1d2c8f'
            )
          ),
          actionButton('btn_persistence_mark_classes',  '', icon = icon('glyphicon glyphicon-plus', lib = 'glyphicon')),
          actionButton('btn_persistence_clear_classes', '', icon = icon("glyphicon glyphicon-fire", lib = 'glyphicon')),
          HTML("&nbsp;&nbsp;&nbsp;&nbsp;"),
          actionButton('btn_persistence_export_image',    '', icon = icon('glyphicon glyphicon-picture',  lib = 'glyphicon')),
          HTML('&nbsp;&nbsp;'),
          # textOutput('log_persistence_available', inline = TRUE),
          # br()
          checkboxInput('switch_persistence_ratio', label = 'Adjust for dissimilar densities', value = FALSE),
          fluidRow(
            h4('Selected homology classes'),
            verbatimTextOutput('log_persistence_selected', placeholder = TRUE)#,
            #tags$head(tags$style("#log_persistence_selected{overflow-y:scroll;}"))
          ),
          fluidRow(
            h4('Marked homology classes'),
            verbatimTextOutput('log_persistence_marked', placeholder = TRUE)#,
            #tags$head(tags$style("#log_persistence_marked{ overflow-y:scroll;}"))
          )
        )
      )
    ),
    ## MIDDLE PANEL: TRAJECTORY DENDROGRAM
    column(
      width = 4,
      style = 'padding-left: 50px',
      # Dendrogram of trajectories clustered by marked homology classes
      tabsetPanel(
        tabPanel(
          title = 'Whole dendrogram',
          fluidRow(
            plotOutput('plot_dendrogram',
                       height = 700,
                       brush  = brushOpts(
                         id        = 'selector_dendrogram',
                         fill      = 'yellow',
                         stroke    = 'yellow',
                         direction = 'y'
                       ))
          ),
          column(
            width = 8,
            actionButton('btn_dendrogram_mark_leaves',                 '', icon = icon('glyphicon glyphicon-plus',       lib = 'glyphicon')),
            actionButton('btn_dendrogram_clear_marked_leaves',         '', icon = icon('glyphicon glyphicon-fire',       lib = 'glyphicon')),
            HTML('&nbsp;&nbsp;&nbsp;&nbsp;'),
            actionButton('btn_dendrogram_zoom',                        '', icon = icon('glyphicon glyphicon-glyphicon glyphicon-zoom-in', lib = 'glyphicon')),
            HTML('&nbsp;&nbsp;&nbsp;&nbsp;'),
            actionButton('btn_trajectories_export_fcs',                '', icon = icon('glyphicon glyphicon-save',       lib = 'glyphicon')),
            actionButton('btn_trajectories_clear_pinned_trajectories', '', icon = icon('glyphicon glyphicon-trash',      lib = 'glyphicon')),
            HTML('&nbsp;&nbsp;&nbsp;&nbsp;'),
            actionButton('btn_dendrogram_export_image',                  '', icon = icon('glyphicon glyphicon-picture',  lib = 'glyphicon'))
          ),
          column(
            width = 4,
            style = 'text-align:  right;
                 padding-right: 5px;
                 margin-top:   15px',
            textOutput('log_pinned_batches_count')
          )
        ),
        tabPanel(
          title = 'Zoom',
          plotOutput('plot_dendrogram_zoom',
                     height = 700,
                     brush  = brushOpts(
                       id        = 'selector_dendrogram_zoom',
                       fill      = 'yellow',
                       stroke    = 'yellow',
                       direction = 'y'
                     )),
          column(
            width = 6,
            actionButton('btn_dendrogram_zoom_mark_leaves',                 '', icon = icon('glyphicon glyphicon-plus',       lib = 'glyphicon')),
            actionButton('btn_dendrogram_zoom_clear_marked_leaves',         '', icon = icon('glyphicon glyphicon-fire',       lib = 'glyphicon')),
          ),
          column(
            width = 6,
            style = 'text-align:  right;
                 padding-right: 5px;
                 margin-top:   15px',
            textOutput('log_pinned_batches_count_zoom')
          )
        )
      ),
      br(), br(), br(),
      # Trajectory group switch (A vs. B)
      radioGroupButtons(
        'btn_trajectories_group',
        label = NULL,
        choiceValues = c('A', 'B'),
        choiceNames  = c('A', 'B'),
        selected = 'A', status = 'default',
        size = 'sm', direction = 'horizontal',
        justified = TRUE, individual = FALSE
      ),
      tags$script("$(\"input:radio[name='btn_trajectories_group'][value='A']\").parent().css('background-color', '#c2dfff');"),
      tags$script("$(\"input:radio[name='btn_trajectories_group'][value='B']\").parent().css('background-color', '#ffb5c9');"),
      sliderInput(
        'slider_dendrogram_leaf_cutoff',
        'Min trajectory count % per leaf',
        min = 0, max = 100,
        value = 10.0,
        step  = 1,
        width = '100%'
      ),
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
          #verbatimTextOutput('log_trajectories_marked.A', placeholder = TRUE)
          htmlOutput('log_trajectories_marked.A', placeholder = TRUE, style = 'padding-right: 5px')
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
          #verbatimTextOutput('log_trajectories_marked.B', placeholder = TRUE)
          htmlOutput('log_trajectories_marked.B', placeholder = TRUE, style = 'padding-left: 5px')
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
        shinycssloaders::withSpinner(
          plotOutput(
            'plot_layout_trajectories',
            height = 500
          ),
          color='#1d2c8f'
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
          actionButton('btn_layout_trajectories_flip_colours', label = '', icon = icon('glyphicon glyphicon-adjust', lib = 'glyphicon')),
          actionButton('btn_layout_trajectories_export_image',           '', icon = icon('glyphicon glyphicon-picture',  lib = 'glyphicon'))
        ),
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
                width = 8,
                style = 'padding-left: 5px',
                selectInput('input_tracked_markers.A', label = 'Group A markers of interest', choices = c(), multiple = TRUE, width = '100%')
              ),
              column(
                width = 4,
                style = 'padding-left:2px;
                     margin-top: 26px;
                    text-align: right',
                actionButton('btn_tracked_markers_remove_trajectories.A', label = '',      icon = icon('glyphicon glyphicon-flash',          lib = 'glyphicon')),
                actionButton('btn_tracked_markers_undo_remove_trajectories.A', label = '', icon = icon('glyphicon glyphicon-step-backward',  lib = 'glyphicon')),
                actionButton('btn_tracked_markers_highlight_segments.A',  label = '',      icon = icon('glyphicon glyphicon-flag',           lib = 'glyphicon')),
                actionButton('btn_tracked_markers_export_image.A',                   '',      icon = icon('glyphicon glyphicon-picture',      lib = 'glyphicon'))
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
                width = 8,
                style = 'padding-left: 5px',
                selectInput('input_tracked_markers.B', label = 'Group B markers of interest', choices = c(), multiple = TRUE, width = '100%')
              ),
              column(
                width = 4,
                style = 'padding-left: 2px;
                   margin-top:  26px;
                   text-align: right',
                actionButton('btn_tracked_markers_remove_trajectories.B', label = '',      icon = icon('glyphicon glyphicon-flash',          lib = 'glyphicon')),
                actionButton('btn_tracked_markers_undo_remove_trajectories.B', label = '', icon = icon('glyphicon glyphicon-step-backward',  lib = 'glyphicon')),
                actionButton('btn_tracked_markers_highlight_segments.B', label = '',       icon = icon('glyphicon glyphicon-flag',           lib = 'glyphicon')),
                actionButton('btn_tracked_markers_export_image.B',                  '',      icon = icon('glyphicon glyphicon-picture',      lib = 'glyphicon'))
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
                width = 9,
                style = 'padding-left: 5px;
                         padding-right: 5px',
                selectInput('input_tracked_populations.A', label = 'Group A populations of interest', choices = c(), multiple = TRUE, width = '100%')
              ),
              column(
                width = 3,
                style = 'padding-left: 2px;
                   margin-top:  26px;
                   text-align: right',
                actionButton('btn_tracked_populations_highlight_segments.A', label = '',       icon = icon('glyphicon glyphicon-flag',           lib = 'glyphicon')),
                actionButton('btn_tracked_populations_export_image.A',  label = '',            icon = icon('glyphicon glyphicon-picture',        lib = 'glyphicon'))
              )
            ),
            fluidRow(
              plotOutput(
                'plot_tracked_populations.A',
                height = 500,
                brush = brushOpts(
                  id     = 'selector_tracked_populations.A',
                  fill   = 'red',
                  stroke = 'red'
                )
              )
            ),
            ## B
            fluidRow(
              column(
                width = 9,
                style = 'padding-left: 5px;
                         padding-right: 5px',
                selectInput('input_tracked_populations.B', label = 'Group B populations of interest', choices = c(), multiple = TRUE, width = '100%')
              ),
              column(
                width = 3,
                style = 'padding-left: 2px;
                   margin-top:  26px;
                   text-align: right',
                actionButton('btn_tracked_populations_highlight_segments.B', label = '',       icon = icon('glyphicon glyphicon-flag',           lib = 'glyphicon')),
                actionButton('btn_tracked_populations_export_image.B',  label = '',            icon = icon('glyphicon glyphicon-picture',        lib = 'glyphicon'))
              )
            ),
            fluidRow(
              
              plotOutput(
                'plot_tracked_populations.B',
                height = 500,
                brush = brushOpts(
                  id     = 'selector_tracked_populations.B',
                  fill   = 'red',
                  stroke = 'red'
                )
              )
              
            )
          )
        )
      )
    )
  )
)