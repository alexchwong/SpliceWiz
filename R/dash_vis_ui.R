ui_vis_diag <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(3,
                selectInput(ns('filterType_diag'), 'Filter Events by', 
                    c("Adjusted P value", "Nominal P value", "Top N results")),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top N results'].",
                        "indexOf(input.filterType_diag) != 0"
                    ),
                    shinyWidgets::sliderTextInput(
                        inputId = ns("pvalT_diag"), 
                        label = "P-value/FDR threshold",
                        choices = c(0.000001, 0.0001, 0.001, 
                            0.01, 0.05, 0.1, 0.2), 
                        selected = 0.05
                    )
                ),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top N results'].",
                        "indexOf(input.filterType_diag) == 0"
                    ),
                    shinyWidgets::sliderTextInput(
                        inputId = ns("topN_diag"), 
                        label = "Top N results",
                        choices = c(10, 20, 50, 100, 200, 300, 
                            500, 1000, 2000, 5000, 10000), 
                        selected = 500
                    )
                ),
                selectInput(ns("EventType_diag"), 
                    "Filter Events by ASE Modality", 
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
                actionButton(ns("output_plot_diag"), 
                    "Generate RStudio ggplot"),
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
                selectInput(ns('filterType_volc'), 'Filter Events by', 
                    c("Adjusted P value", "Nominal P value", "Top N results")),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top N results'].",
                        "indexOf(input.filterType_volc) != 0"
                    ),
                    shinyWidgets::sliderTextInput(
                        inputId = ns("pvalT_volc"), 
                        label = "P-value/FDR threshold",
                        choices = c(0.000001, 0.0001, 0.001, 
                            0.01, 0.05, 0.1, 0.2), 
                        selected = 0.05
                    )
                ),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top N results'].",
                        "indexOf(input.filterType_volc) == 0"
                    ),
                    shinyWidgets::sliderTextInput(
                        inputId = ns("topN_volc"), 
                        label = "Top N results",
                        choices = c(10, 20, 50, 100, 200, 300, 
                            500, 1000, 2000, 5000, 10000), 
                        selected = 500
                    )
                ),
                selectInput(ns("EventType_volc"), 
                    "Filter Events by ASE Modality", 
                    width = '100%', multiple = TRUE,
                    choices = c("IR", "MXE", "SE", "AFE", "ALE", 
                        "A5SS", "A3SS")),
                shinyWidgets::switchInput(ns("facet_volc"), 
                    label = "Facet by ASE Modality", 
                    labelWidth = "150px"),
                shinyWidgets::switchInput(ns("adjP_volc"), 
                    label = "Use Adjusted P values", 
                    value = TRUE, labelWidth = "100px"),
                shinyWidgets::switchInput(ns("NMD_volc"), 
                    label = "NMD Mode", labelWidth = "100px"),
                # shinySaveButton(ns("saveplot_volc"), 
                    # "Save Plot as PDF", "Save Plot as PDF...", 
                    # filetype = list(PDF = "pdf")),
                actionButton(ns("output_plot_volc"), 
                    "Generate RStudio ggplot"),
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
                selectInput(ns('filterType_heat'), 'Filter Events by', 
                    c("Adjusted P value", "Nominal P value", "Top N results")),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top N results'].",
                        "indexOf(input.filterType_heat) != 0"
                    ),
                    shinyWidgets::sliderTextInput(
                        inputId = ns("pvalT_heat"), 
                        label = "P-value/FDR threshold",
                        choices = c(0.000001, 0.0001, 0.001, 
                            0.01, 0.05, 0.1, 0.2), 
                        selected = 0.05
                    )
                ),
                conditionalPanel(ns = ns,
                    condition = paste0(
                        "['Top N results'].",
                        "indexOf(input.filterType_heat) == 0"
                    ),
                    shinyWidgets::sliderTextInput(
                        inputId = ns("topN_heat"), 
                        label = "Top N results",
                        choices = c(10, 20, 50, 100, 200, 300, 
                            500, 1000, 2000, 5000, 10000), 
                        selected = 500
                    )
                ),
                selectInput(ns('secondFilter_heat'), 
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
                        "indexOf(input.secondFilter_heat) == 0"
                    ),
                    selectInput(ns('GO_heat'), 'Filter by GO category', 
                        c("(none)")),
                ),
                selectInput(ns("anno_col_heat"), 
                    "Display Annotation Categories", 
                    width = '100%', multiple = TRUE,
                    choices = c()),
                selectInput(ns("anno_col_heat_sort"), "Sort by Category", 
                    width = '100%', choices = "(none)"),
                shinyWidgets::radioGroupButtons(ns("anno_col_heat_sort_order"), 
                    label = "Sort order", justified = FALSE,
                    choiceNames = c("Ascending", "Descending"), 
                    choiceValues = c(FALSE, TRUE),
                    checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                ),
                shinyWidgets::sliderTextInput(ns("slider_num_events_heat"), 
                    "Number of Top Events", 
                    choices = c(5, 10,25,50,100,200,500), 
                    selected = 25),
                shinyWidgets::radioGroupButtons(ns("mode_heat"), 
                    label = "Mode", justified = FALSE,
                    choices = c("PSI", "Logit", "Z-score"), 
                    checkIcon = list(yes = icon("ok", lib = "glyphicon"))
                ),                   
                selectInput(ns('color_heat'), 'Palette', 
                    c("RdBu", "BrBG", "PiYG", "PRGn", "PuOr", 
                        "RdGy", "RdYlBu", "RdYlGn", "Spectral"),
                    selected = "RdYlBu"
                ),
                # shinySaveButton(ns("saveplot_heat"), 
                    # "Save Plot as PDF", "Save Plot as PDF...", 
                    # filetype = list(PDF = "pdf"))
                # actionButton(ns("output_plot_heat"), 
                    # "Generate RStudio ggplot"),
            ),
            column(9, 
                textOutput(ns("warning_heat")),
                plotlyOutput(ns("plot_heat"), height = "800px"),
            )
        )    
    )
}