shiny_ui <- fluidPage(
  tags$style(HTML("#app_title{font-family: monospace; font-weight: bold}")),
  tags$style(type = "text/css", ".irs-slider {width: 15px; height: 30px; top: 22px;};"),
  
  titlePanel(
    h3(id = "app_title", "tviblindi"),
    windowTitle = "tviblindi"
  ),
  br(),
      mainPanel(
          style = "width:100%",
          column(
              width = 4,
              
              tabsetPanel(
                tabPanel(
                  "Terminal nodes",
                  fluidRow(
                    plotOutput("term_plot", height = 500,
                               brush = brushOpts(
                                 id     = "term_brush",
                                 fill   = "yellow",
                                 stroke = "yello"
                               )
                    )
                  ),
                  
                  ## Button: mark selected terminal nodes
                  actionButton("term_btn_add", "", icon = icon("glyphicon glyphicon-plus", lib = "glyphicon")),
                  
                  ## Button: clear marked teminal nodes
                  actionButton("term_btn_clear", "", icon = icon("glyphicon glyphicon-fire", lib = "glyphicon")),
                  
                  HTML("&nbsp;&nbsp;&nbsp;&nbsp;"),
                  
                  ## Button: update walks by selected terminal nodes
                  actionButton("term_btn_update", "", icon = icon(" glyphicon glyphicon-thumbs-up", lib = "glyphicon")),
                  
                  ## Selected terminal nodes log
                  fluidRow(
                    h4("Selected terminal nodes"),
                    verbatimTextOutput("term_brush_info", placeholder = TRUE)
                  ),
                  fluidRow(
                    h4("Terminal nodes marked for clustering"),
                    verbatimTextOutput("term_marked_info", placeholder = TRUE)
                  )
                ),
                tabPanel(
                  "Persistence",
                  ## Interactive persistence diagram
                  fluidRow(
                    plotOutput("pers_plot", height = 500,
                               brush = brushOpts(
                                 id     = "pers_brush",
                                 fill   = "yellow",
                                 stroke = "yellow"
                               )
                    )
                  ),
                  
                  ## Button: mark selected persistence diagram points
                  actionButton("pers_btn_add", "", icon = icon("glyphicon glyphicon-plus", lib = "glyphicon")),
                  
                  ## Button: clear marked persistence diagram points
                  actionButton("pers_btn_clear", "", icon = icon("glyphicon glyphicon-fire", lib = "glyphicon")),
                  
                  HTML("&nbsp;&nbsp;"),
                  
                  textOutput("walks_status", inline = TRUE),
                  
                  ## Selected homology classes log
                  fluidRow(
                    h4("Selected points"),
                    verbatimTextOutput("pers_brush_info", placeholder = TRUE)
                  ),
                  
                  ## Marked homology classes log
                  fluidRow(
                    h4("Marked points"),
                    verbatimTextOutput("pers_marked_info", placeholder = TRUE)
                  )    
                )
              )
              
          ),
          
          column(
              width = 4,
              style = "padding-left:50px;",
              
              ## Interactive pathways dendrogram
              fluidRow(
                  plotOutput("dendro_plot", height = 800,
                             brush = brushOpts(
                                 id        = "dendro_brush",
                                 fill      = "yellow",
                                 stroke    = "yellow",
                                 direction = "y"
                             ))
              ),
              
              ## Pathway percentage threshold selector
              sliderInput("dendro_perc", "Min pathway count % per leaf",
                          min = 0, max = 100,
                          value = 10.0, step = 1),
              
             ## Button switch: pathway category (A vs. B)
             radioGroupButtons("dendro_btn_categ", label = NULL, choiceValues = c("A", "B"), choiceNames = c("A", "B"),
                               selected = "A", status = "default", size = "sm",
                               direction = "horizontal", justified = TRUE, individual = FALSE),
             
              column(width = 6,
                     ## Button: mark selected pathways
                     actionButton("dendro_btn_add", "", icon = icon("glyphicon glyphicon-plus", lib = "glyphicon")),
                     
                     ## Button: clear marked pathways
                     actionButton("dendro_btn_clear", "", icon = icon("glyphicon glyphicon-fire", lib = "glyphicon")),
                     
                     HTML("&nbsp;&nbsp;&nbsp;&nbsp;"),
                     
                     ## Button: open dialog for saving output .fcs file
                     actionButton("dendro_btn_save", "", icon = icon("glyphicon glyphicon-floppy-disk", lib = "glyphicon")),     
                     
              ),
              column(width = 6,
                     style = "text-align:right; padding-right:5px; margin-top: 15px;",
                     textOutput("save_count_info")
              ),

             br(),
             br(),
             
              ## Selected pathways log
              fluidRow(
                  width = 12,
                  h4("Selected dendrogram nodes"),
                  verbatimTextOutput("dendro_brush_info", placeholder = TRUE)
              ),
              
              fluidRow(
                style = "padding-left:0px; padding-right:0px",
                ## Marked pathways in category A log
                column(
                  width = 6,
                  style = "padding-left:0px; padding-right:5px",
                  
                  h4("Marked pathways in A"),
                  
                  column(width = 6,
                         ## Button: pin marked pathways for addition to output .fcs file
                         actionButton("dendro_btn_append.A", "", icon = icon("glyphicon glyphicon-pushpin", lib = "glyphicon")),
                  ),
                  
                  column(width = 6,
                         textOutput("marked_A_deleted_info")
                  ),
                  
                  br(), br(),
                  
                  textOutput("marked_A_info")
                ),
                
                column(
                  width = 6,
                  style = "padding-left:0px; padding-right:5px",
                  
                  h4("Marked pathways in B"),
                  
                  column(width = 6,
                         ## Button: pin marked pathways for addition to output .fcs file
                         actionButton("dendro_btn_append.B", "", icon = icon("glyphicon glyphicon-pushpin", lib = "glyphicon")),
                  ),
                  
                  column(width = 6,
                         textOutput("marked_B_deleted_info")
                  ),
                  
                  br(), br(),
                  
                  textOutput("marked_B_info")
                )
              ),
              
              br(), br(), br()
          ),
          column(
              width = 4,
              style = "padding-left:50px;",
              
              ## Data layout projection
              fluidRow(
                  plotOutput("layout_plot", height = 500)
              ),
              
              ## Checkbox: lazy plotting
              fluidRow(
                column(width = 6,
                       checkboxInput("layout_lazy", "Lazy plotting", value = TRUE, width = NULL)),
                column(width = 6,
                       style = "text-align:right; padding-right:5px; margin-top: 15px;",
                       actionButton("layout_btn_flip_colours", label = "", icon = icon("glyphicon glyphicon-adjust", lib = "glyphicon")))
              ),
              
              fluidRow(
                  column(width = 6,
                         style = "padding-left:5px;",
                         ## Marker selector
                         selectInput("marker_selector.A", label = "Group A markers of interest", choices = c("foo"), multiple = TRUE)  
                  ),
                  column(width = 2,
                         style = "padding-left:2px;",
                         ## Scaling exponent
                         numericInput("scaling_exponent.A", label = "Scale exp", min = 0.1, max = 10, value = 1, step = 0.1)
                  ),
                  column(width = 2,
                         style = "padding-left:2px;",
                         ## Segment number selector
                         numericInput("n_segments.A", label = "Segs", min = 3, max = 1000, value = 10)
                  ),
                  column(width = 2,
                         style = "padding-left:2px; margin-top: 26px;",
                         ## Pathway zapper
                         actionButton("expression_zap.A", label = "", icon = icon("glyphicon glyphicon-flash", lib = "glyphicon"))
                  )
              ),
              
              ## Marker expression plot A
              fluidRow(
                  plotOutput("expression_plot.A", height = 500,
                             brush = brushOpts(
                               id        = "expression_brush.A",
                               fill      = "red",
                               stroke    = "red"
                             ))
              ),
              fluidRow(
                  column(width = 6,
                         style = "padding-left:5px;",
                         ## Marker selector
                         selectInput("marker_selector.B", label = "Group B markers of interest", choices = c("foo"), multiple = TRUE)  
                  ),
                  column(width = 2,
                         style = "padding-left:2px;",
                         ## Scaling exponent
                         numericInput("scaling_exponent.B", label = "Scale exp", min = 0.1, max = 10, value = 1, step = 0.1)
                  ),
                  column(width = 2,
                         style = "padding-left:2px;",
                         ## Segment number selector
                         numericInput("n_segments.B", label = "Segs", min = 3, max = 1000, value = 10)
                  ),
                  column(width = 2,
                         style = "padding-left:2px; margin-top: 26px;",
                         ## Pathway zapper
                         actionButton("expression_zap.B", label = "", icon = icon("glyphicon glyphicon-flash", lib = "glyphicon"))
                  )
              ),
              ## Marker expression plot B
              fluidRow(
                  plotOutput("expression_plot.B", height = 500,
                             brush = brushOpts(
                               id        = "expression_brush.B",
                               fill      = "red",
                               stroke    = "red"
                             ))
              )
          )
      )
  )
