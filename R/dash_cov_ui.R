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
            column(4, 
                wellPanel(style = "overflow-y:scroll; max-height: 800px",
                    selectInput(ns('event_norm_cov'), 'Normalize by Event', 
                        width = '100%', choices = c("(none)")),
                    selectInput(ns('condition_cov'), 'Condition variable:', 
                        width = '100%', choices = c("(Individual Samples)")),
                    wellPanel(style = "overflow-y:scroll; max-height: 300px",
                        rHandsontableOutput(ns("track_table"))
                    ),
                    
                    conditionalPanel(ns = ns,
                        condition = "['(Individual Samples)'].indexOf(input.condition_cov) < 0",
                        # selectInput(ns('diff_stat'), 
                            # 'Differential coverage track', 
                            # width = '100%', choices = c("t-test")),                        
                        shinyWidgets::prettyRadioButtons(
                            inputId = ns("diff_stat"),
                            label = "Measure difference by:", 
                            choiceNames = c("t-test", "none"),
                            choiceValues = c("t-test", "none"),
                            inline = TRUE, 
                            status = "danger",
                            fill = TRUE
                        ),
                        selectInput(ns('diffA'), 'Contrasting category A', 
                            width = '100%', choices = c("(none)")),
                        selectInput(ns('diffB'), 'Contrasting category B', 
                            width = '100%', choices = c("(none)")),
                    ),
                    shinyWidgets::sliderTextInput(ns("slider_num_plotRes"), 
                        "Plot resolution", 
                        choices = c(250, 500, 1000, 2000, 5000, 10000),
                        selected = 1000),
                    shinyWidgets::prettyRadioButtons(
                        inputId = ns("plot_ribbon"),
                        label = "Variance by:", 
                        choiceNames = c("SD", "SEM", "95% CI", "none"),
                        choiceValues = c("sd", "sem", "ci", "none"),
                        inline = TRUE, 
                        status = "danger",
                        fill = TRUE
                    ),
                    shinyWidgets::materialSwitch(
                       inputId = ns("plot_Jn_cov"),
                       label = "Plot Junctions", right = TRUE,
                       value = FALSE, status = "success"
                    ),
                    shinyWidgets::materialSwitch(
                       inputId = ns("normalizeCov"),
                       label = "Normalize Coverage", right = TRUE,
                       value = FALSE, status = "success"
                    ),
                    shinyWidgets::materialSwitch(
                       inputId = ns("plot_key_iso"),
                       label = "Display Event Isoforms only", right = TRUE,
                       value = FALSE, status = "success"
                    ),
                    shinyWidgets::materialSwitch(
                       inputId = ns("condense_cov"),
                       label = "Condense Annotation by Gene", right = TRUE,
                       value = FALSE, status = "success"
                    ),
                    shinyWidgets::materialSwitch(
                       inputId = ns("selTr_cov"),
                       label = "Select transcripts", right = TRUE,
                       value = FALSE, status = "success"
                    ),
                    shinyWidgets::materialSwitch(
                       inputId = ns("exonMode_cov"),
                       label = "Exon Plot Mode", right = TRUE,
                       value = FALSE, status = "success"
                    )
                ),
                conditionalPanel(ns = ns,
                    condition = "input.selTr_cov == true",                
                    wellPanel(style = "overflow-y:scroll; max-height: 250px",
                        rHandsontableOutput(ns("transcripts_lookup"))
                    )
                ),
                conditionalPanel(ns = ns,
                    condition = "input.exonMode_cov == true",                
                    wellPanel(style = "overflow-y:scroll; max-height: 250px",
                        rHandsontableOutput(ns("exons_lookup"))
                    )
                ),
                vis_ggplot_UI(ns("covSave"), "Save interactive plot as PDF"),
                br(),
                vis_ggplot_UI(ns("covExonSave"), "Save exon-centric plot as PDF"),
            ),
            column(8, 
                plotlyOutput(ns("plot_cov"), height = "600px"),
                conditionalPanel(ns = ns,
                    condition = "input.exonMode_cov == true",
                    plotOutput(ns("stillplot_cov"), height = "600px")
                )
            )
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