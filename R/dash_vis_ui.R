vis_ggplot_UI <- function(id, label = "Export as ggplot") {
    ns <- NS(id)
    ui_toggle_wellPanel_modular(
        inputId = "ddb_ggsave",
        id = id,
        title = "Save Plot to PDF",
        color = "default",
        icon = icon("folder-open", lib = "font-awesome"), 
        wellPanel(
            selectInput(ns('units'), 'Units:', 
                choices = c("in", "cm", "mm", "px"),
                selected = "in"
            ),
            numericInput(ns("wt"), "Width:", 8),
            numericInput(ns("ht"), "Height:", 6),
            shinySaveButton(ns("saveplot"),
                "Save", "Save Plot as PDF...", 
                buttonType = "primary",
                filetype = list(PDF = "pdf")
            ),
        )
    )
}

visFilter_UI <- function(id, label = "Filter events by") {
    ns <- NS(id)
    wellPanel(
        h5(label),
        selectInput(ns('vF_filterType'), 'Filter Events by', 
            choices = c(
                "Adjusted P value", 
                "Nominal P value", 
                "Top events by p-value",
                "Highlighted events"
            ),
            selected = "Adjusted P value"
        ),
        conditionalPanel(ns = ns,
            condition = paste0(
                "['Top events by p-value'].",
                "indexOf(input.vF_filterType) != 0"
            ),
            shinyWidgets::sliderTextInput(
                inputId = ns("vF_pvalT"), 
                label = "P-value/FDR threshold",
                choices = c(0.000001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 1), 
                selected = 0.05
            )
        ),
        conditionalPanel(ns = ns,
            condition = paste0(
                "['Top events by p-value'].",
                "indexOf(input.vF_filterType) == 0"
            ),
            shinyWidgets::sliderTextInput(
                inputId = ns("vF_topN"), 
                label = "Number of top events",
                choices = c(10, 20, 50, 100, 200, 300, 
                    500, 1000, 2000, 5000, 10000), 
                selected = 500
            )
        ),
        selectInput(ns("vF_EventType"), 
            "Filter Events by ASE Modality", 
            width = '100%', multiple = TRUE,
            choices = c("IR", "MXE", "SE", "AFE", "ALE", 
                "A5SS", "A3SS")), br(),
        actionButton(ns("vF_reset"), "Reset to default"),
    )
}
        

ui_vis_diag <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(3,
                visFilter_UI(ns("scatter")),
                selectInput(ns('variable_diag'), 'Variable', 
                    c("(none)")),
                selectInput(ns('nom_diag'), 'X-axis condition', 
                    c("(none)")),
                selectInput(ns('denom_diag'), 'Y-axis condition', 
                    c("(none)")),
                shinyWidgets::materialSwitch(
                   inputId = ns("NMD_diag"),
                   label = "NMD Mode", right = TRUE,
                   value = FALSE, status = "success"
                ),
                vis_ggplot_UI(ns("scatterSave")),
                # actionButton(ns("output_plot_diag"), 
                    # "Generate RStudio ggplot"),
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
                visFilter_UI(ns("volcano")),
                shinyWidgets::materialSwitch(
                   inputId = ns("facet_volc"),
                   label = "Facet by ASE Modality", right = TRUE,
                   value = FALSE, status = "success"
                ),            shinyWidgets::materialSwitch(
                   inputId = ns("adjP_volc"),
                   label = "Plot adjusted P values", right = TRUE,
                   value = TRUE, status = "success"
                ),            shinyWidgets::materialSwitch(
                   inputId = ns("NMD_volc"),
                   label = "NMD Mode", right = TRUE,
                   value = FALSE, status = "success"
                ),
                vis_ggplot_UI(ns("volcanoSave")),
                # shinySaveButton(ns("saveplot_volc"), 
                    # "Save Plot as PDF", "Save Plot as PDF...", 
                    # filetype = list(PDF = "pdf")),
                # actionButton(ns("output_plot_volc"), 
                    # "Generate RStudio ggplot"),
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
                visFilter_UI(ns("heatmap")),
                selectInput(ns('GO_heat'), 'Filter by GO category', 
                        c("(none)")),
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
                    "Maximum number of rows", 
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
                vis_ggplot_UI(ns("heatSave")),
                # shinySaveButton(ns("saveplot_heat"), 
                    # "Save Plot as PDF", "Save Plot as PDF...", 
                    # filetype = list(PDF = "pdf"))
                actionButton(ns("output_plot_heat"), 
                    "Generate RStudio ggplot"),
            ),
            column(9, 
                textOutput(ns("warning_heat")),
                plotlyOutput(ns("plot_heat"), height = "800px"),
            )
        )    
    )
}