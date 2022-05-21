ui_vis_diag <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(3,
                shinyWidgets::sliderTextInput(
                    inputId = ns("number_events_diag"), 
                    label = "Number of Top Events",
                    choices = c(100, 200, 500, 1000, 2000 , 5000,
                        10000, 20000, 50000, 100000, 200000, 500000), 
                    selected = 1000),
                selectInput(ns("EventType_diag"), "Splice Type", 
                    width = '100%', multiple = TRUE,
                    choices = c("IR", "MXE", "SE", "AFE", "ALE", 
                        "A5SS", "A3SS")),
                selectInput(ns('variable_diag'), 'Variable', 
                    c("(none)")),
                selectInput(ns('nom_diag'), 'X-axis condition', 
                    c("(none)")),
                selectInput(ns('denom_diag'), 'Y-axis condition', 
                    c("(none)")),
                shinyWidgets::switchInput(ns("NMD_diag"), 
                    label = "NMD Mode", labelWidth = "100px"),
                # shinySaveButton(ns("saveplot_diag"), 
                    # "Save Plot as PDF", "Save Plot as PDF...", 
                    # filetype = list(PDF = "pdf")),
                actionButton(ns("clear_diag"), "Clear settings"), br(), br(),
                textOutput(ns("warning_diag"))
            ),
            column(9,
                plotlyOutput(ns("plot_diag"), height = "800px")
            )
        )
    )
}

ui_vis_volcano <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(3,    
                shinyWidgets::sliderTextInput(
                    inputId = ns("number_events_volc"), 
                    label = "Number of Top Events",
                    choices = c(100, 200, 500, 1000, 2000, 5000, 10000,
                        20000, 50000, 100000, 200000, 500000),
                    selected = 1000),
                selectInput(ns("EventType_volc"), "Splice Type", 
                    width = '100%', multiple = TRUE,
                    choices = c("IR", "MXE", "SE", "AFE", "ALE", 
                        "A5SS", "A3SS")),
                shinyWidgets::switchInput(ns("facet_volc"), 
                    label = "Facet by Type", labelWidth = "150px"),
                shinyWidgets::switchInput(ns("adjP_volc"), 
                    label = "Multiple Testing", 
                    value = TRUE, labelWidth = "100px"),
                shinyWidgets::switchInput(ns("NMD_volc"), 
                    label = "NMD Mode", labelWidth = "100px"),
                # shinySaveButton(ns("saveplot_volc"), 
                    # "Save Plot as PDF", "Save Plot as PDF...", 
                    # filetype = list(PDF = "pdf")),
                actionButton(ns("clear_volc"), "Clear settings"), br(), br(),
                textOutput(ns("warning_volc"))
            ),
            column(9,
                plotlyOutput(ns("plot_volc"), height = "800px")
            )
        )
    )
}

ui_vis_heatmap <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(3, 
                shinyWidgets::radioGroupButtons(ns("select_events_heat"), 
                    label = 
                        "Select Events from Differential Expression Results",
                    justified = FALSE,
                    choices = c("Highlighted", "Top N Filtered Results", 
                        "Top N All Results"), 
                    checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                ),
                selectInput(ns("anno_col_heat"), "Sample Annotation Rows", 
                    width = '100%', multiple = TRUE,
                    choices = c()),
                shinyWidgets::sliderTextInput(ns("slider_num_events_heat"), 
                    "Num Events", choices = c(5, 10,25,50,100,200,500), 
                    selected = 25),
                shinyWidgets::radioGroupButtons(ns("mode_heat"), 
                    label = "Mode", justified = FALSE,
                    choices = c("PSI", "Logit", "Z-score"), 
                    checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                ),
                # shinyWidgets::switchInput("NMD_heat", label = "NMD Mode", 
                    # labelWidth = "100px"),                    
                selectInput(ns('color_heat'), 'Palette', 
                    c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", 
                        "RdGy", "RdYlBu", "RdYlGn", "Spectral")
                ),
                # shinySaveButton(ns("saveplot_heat"), 
                    # "Save Plot as PDF", "Save Plot as PDF...", 
                    # filetype = list(PDF = "pdf"))                
            ),
            column(9, 
                textOutput(ns("warning_heat")),
                plotlyOutput(ns("plot_heat"), height = "800px"),
            )
        )    
    )
}