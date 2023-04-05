ui_cov <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(#style='height:25vh',
            column(4, 
                textOutput(ns("warning_cov")),
                div(style=.cov_ui_inline(250), 
                    selectizeInput(ns('genes_cov'), 'Genes', 
                        choices = "(none)")
                ),
                div(style=.cov_ui_inline(350), 
                    selectizeInput(ns('events_cov'), 'Events', 
                        choices = "(none)")
                ), br(),
                actionButton(ns("refresh_coverage"), "Refresh Plot")
            ),
            column(4,
                div(style=.cov_ui_inline(80),
                    selectInput(ns("chr_cov"), label = "Chr", 
                        choices = c("(none)"), selected = "(none)")),
                div(style=.cov_ui_inline(120),
                    numericInput(ns("start_cov"), label = "Left", 
                        value = c(""))),
                div(style=.cov_ui_inline(120),
                    numericInput(ns("end_cov"), label = "Right", 
                        value = c(""))),           
                br(),
                shinyWidgets::actionBttn(ns("zoom_out_cov"), 
                    style = "material-circle", 
                    color = "danger",icon = icon("minus")),
                div(style=.cov_ui_inline(50, "center", padding = 25),
                    textOutput(ns("label_zoom_cov"))),
                shinyWidgets::actionBttn(ns("zoom_in_cov"), 
                    style = "material-circle", 
                    color = "danger", icon = icon("plus")),
                div(style=.cov_ui_inline(alignment = "center", padding = 15),
                    shinyWidgets::radioGroupButtons(ns("strand_cov"), 
                        label = "RNA-seq Strand Coverage", justified = FALSE,
                        choices = c("*", "+", "-"), 
                        checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                    )
                )#, br(),
                # shinyWidgets::radioGroupButtons(ns("graph_mode_cov"), 
                    # label = "Graph Mode", justified = FALSE,
                    # choices = c("Pan", "Zoom", "Movable Labels"), 
                    # checkIcon = list(yes = icon("ok", lib = "glyphicon"))),
            ),
            column(4,
                selectInput(ns('modeFilter_COV'), 
                    'Display events from:', 
                    c(
                        "All filtered events", 
                        "Highlighted (selected) events", 
                        "Top Gene Ontology Categories"
                    )
                ),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top Gene Ontology Categories'].",
                        "indexOf(input.modeFilter_COV) == 0"
                    ),
                    selectInput(ns('GOterm_COV'), 'Filter by GO category', 
                        c("(none)")),
                ),
                shinyWidgets::sliderTextInput(ns("slider_num_events_cov"), 
                    "Number of Top Events", 
                    choices = c(5, 10, 25, 50, 100, 200, 500),
                    selected = 25) 
            )
        ),
        fluidRow(
            column(3, 
                wellPanel(style = "overflow-y:scroll; max-height: 800px",
                    selectInput(ns('mode_cov'), 'View', width = '100%',
                        choices = c("Individual", "By Condition")),
                    selectInput(ns('event_norm_cov'), 'Normalize Event', 
                        width = '100%', choices = c("(none)")),
                    conditionalPanel(ns = ns,
                        condition = "['By Condition'].indexOf(input.mode_cov) >= 0",
                        selectInput(ns('condition_cov'), 'Condition', 
                            width = '100%', choices = c("(none)"))
                    ),
                    selectInput(ns('track1_cov'), 'Track 1', width = '100%',
                        choices = c("(none)")),
                    selectInput(ns('track2_cov'), 'Track 2', width = '100%',
                        choices = c("(none)")),
                    selectInput(ns('track3_cov'), 'Track 3', width = '100%',
                        choices = c("(none)")),
                    selectInput(ns('track4_cov'), 'Track 4', width = '100%',
                        choices = c("(none)")),
                    shinyWidgets::switchInput(ns("plot_Jn_cov"), 
                        label = "Plot Junctions", labelWidth = "150px"),
                    shinyWidgets::switchInput(ns("stack_tracks_cov"), 
                        label = "Stack Traces", labelWidth = "150px"),
                    shinyWidgets::switchInput(ns("pairwise_t_cov"), 
                        label = "Pairwise t-test", labelWidth = "150px"),
                    shinyWidgets::switchInput(ns("plot_key_iso"), 
                        label = "Display Key Isoforms", labelWidth = "150px"),
                    shinyWidgets::switchInput(ns("condense_cov"), 
                        label = "Condensed Tracks", labelWidth = "150px"),
                    # shinySaveButton(ns("saveplot_cov"), 
                        # "Save Plot as PDF", "Save Plot as PDF...", 
                        # filetype = list(PDF = "pdf")),
                )
            ),
            column(9, plotlyOutput(ns("plot_cov"), height = "800px"))
        )    
    )
}

ui_cov_new <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(4, 
                textOutput(ns("warning_cov")),
                div(style=.cov_ui_inline(250), 
                    selectizeInput(ns("genes_cov"), 'Locate by Gene', 
                        choices = "(none)")
                ),
                div(style=.cov_ui_inline(350), 
                    selectizeInput(ns("events_cov"), 'Locate by Event', 
                        choices = "(none)")
                ) #, br(),
                # actionButton(ns("refresh_coverage"), "Refresh Plot")
            ),
            column(4,
                div(style=.cov_ui_inline(80),
                    selectInput(ns("chr_cov"), label = "Chr", 
                        choices = c("(none)"), selected = "(none)")),
                div(style=.cov_ui_inline(120),
                    numericInput(ns("start_cov"), label = "Left", 
                        value = c(""))),
                div(style=.cov_ui_inline(120),
                    numericInput(ns("end_cov"), label = "Right", 
                        value = c(""))),           
                br(),
                shinyWidgets::actionBttn(ns("zoom_out_cov"), 
                    style = "material-circle", 
                    color = "danger",icon = icon("minus")),
                div(style=.cov_ui_inline(50, "center", padding = 25),
                    textOutput(ns("label_zoom_cov"))),
                shinyWidgets::actionBttn(ns("zoom_in_cov"), 
                    style = "material-circle", 
                    color = "danger", icon = icon("plus")),
                div(style=.cov_ui_inline(alignment = "center", padding = 15),
                    shinyWidgets::radioGroupButtons(ns("strand_cov"), 
                        label = "RNA-seq Strand Coverage", justified = FALSE,
                        choices = c("*", "+", "-"), 
                        checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                    )
                )#, br(),
                # shinyWidgets::radioGroupButtons(ns("graph_mode_cov"), 
                    # label = "Graph Mode", justified = FALSE,
                    # choices = c("Pan", "Zoom", "Movable Labels"), 
                    # checkIcon = list(yes = icon("ok", lib = "glyphicon"))),
            ),
            column(4,
                selectInput(ns('modeFilter_COV'), 
                    'Display available events by:', 
                    c(
                        "All filtered events", 
                        "Highlighted (selected) events", 
                        "Top Gene Ontology Categories"
                    )
                ),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top Gene Ontology Categories'].",
                        "indexOf(input.modeFilter_COV) == 0"
                    ),
                    selectInput(ns('GOterm_COV'), 'Filter by GO category', 
                        c("(none)")),
                ),
                shinyWidgets::sliderTextInput(ns("slider_num_events_cov"), 
                    "Number of Top Events", 
                    choices = c(5, 10, 25, 50, 100, 200, 500),
                    selected = 25) 
            )
        ),
        fluidRow(
            column(3, 
                wellPanel(style = "overflow-y:scroll; max-height: 800px",
                    selectInput(ns('event_norm_cov'), 'Normalize by Event', 
                        width = '100%', choices = c("(none)")),
                    # selectInput(ns('mode_cov'), 'Tracks by:', width = '100%',
                        # choices = c("Individual samples", "By Condition")),
                    # conditionalPanel(ns = ns,
                        # condition = "['By Condition'].indexOf(input.mode_cov) >= 0",
                        # selectInput(ns('condition_cov'), 'Condition variable:', 
                            # width = '100%', choices = c("(none)"))
                    # ),
                    selectInput(ns('condition_cov'), 'Condition variable:', 
                        width = '100%', choices = c("(Individual Samples)")),
                    wellPanel(style = "overflow-y:scroll; max-height: 300px"
                        rHandsontableOutput(ns("track_table"))
                    ),
                    
                    conditionalPanel(ns = ns,
                        condition = "['By Condition'].indexOf(input.mode_cov) >= 0",
                        selectInput(ns('diff_stat'), 
                            'Differential coverage track', 
                            width = '100%', choices = c("t-test")),                        
                        selectInput(ns('diffA'), 'Contrasting category A', 
                            width = '100%', choices = c("(none)")),
                        selectInput(ns('diffB'), 'Contrasting category B', 
                            width = '100%', choices = c("(none)")),
                    ),

                    shinyWidgets::materialSwitch(
                       inputId = ns("plot_Jn_cov"),
                       label = "Plot Junctions", right = TRUE,
                       value = FALSE, status = "success"
                    )
                    shinyWidgets::materialSwitch(
                       inputId = ns("plot_key_iso"),
                       label = "Display Event Isoforms only", right = TRUE,
                       value = FALSE, status = "success"
                    )
                    shinyWidgets::materialSwitch(
                       inputId = ns("condense_cov"),
                       label = "Condense Annotation by Gene", right = TRUE,
                       value = FALSE, status = "success"
                    )
                )
            ),
            column(9, plotlyOutput(ns("plot_cov"), height = "800px"))
        )    
    )
}

.cov_ui_inline <- function(width, alignment = "top", padding) {
    string <- paste0(
        "display: inline-block;",
        "vertical-align:", alignment, ";"
    )
    if(!missing(width)) {
        string <- paste0(string,
            "width: ", as.character(width), "px;"
        )
    }
    if(!missing(padding)) {
        string <- paste0(string,
            "padding:", as.character(padding), "px"
        )
    }
    return(string)
}