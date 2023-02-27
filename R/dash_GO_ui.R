ui_GO <- function(id) {
    ns <- NS(id)
    fluidRow(
        .ui_notice(),
        column(3,
            textOutput(ns("warning_GO")),
            selectInput(ns('category_GO'), 'GO Category', 
                c("Biological Pathway", "Molecular Function", 
                    "Cellular Compartment")),
            selectInput(ns('threshType_GO'), 'Filter Events by', 
                c("Adjusted P value", "Nominal P value", 
                    "Highlighted events", "Top events by p-value"
                )),
            conditionalPanel(ns = ns,
                condition = "['Top events by p-value'].indexOf(input.threshType_GO) != 0",
                shinyWidgets::sliderTextInput(
                    inputId = ns("pvalT_GO"), 
                    label = "P-value/FDR threshold",
                    choices = c(0.000001, 0.0001, 0.001, 
                        0.01, 0.05, 0.1, 0.2), 
                    selected = 0.05
                )
            ),
            conditionalPanel(ns = ns,
                condition = "['Top events by p-value'].indexOf(input.threshType_GO) == 0",
                shinyWidgets::sliderTextInput(
                    inputId = ns("topN_GO"), 
                    label = "Number of top events",
                    choices = c(10, 20, 50, 100, 200, 300, 
                        500, 1000, 2000, 5000, 10000), 
                    selected = 500
                )
            ),
            selectInput(ns("EventType_GO"), "Filter Events by ASE Modality", 
                width = '100%', multiple = TRUE,
                choices = c("IR", "MXE", "SE", "AFE", "ALE", 
                    "A5SS", "A3SS")
            ),
            shinyWidgets::radioGroupButtons(ns("direction_GO"), 
                label = "Regulation direction", justified = FALSE,
                choices = c("Up", "Down", "Both"), 
                checkIcon = list(yes = icon("ok", lib = "glyphicon"))
            ),
            selectInput(ns('universe_GO'), 'Background Genes based on', 
                c("All ASE events", "Selected ASE Modality", "All Genes")),
            actionButton(ns("perform_GO"), "Perform Gene Ontology Analysis"),
            br(), br(),
            shinySaveButton(ns("GO_export_geneId"), 
                "Save Gene ID's to file", "Save Gene ID's to file...", 
                filetype = list(txt = "txt")),
            shinySaveButton(ns("GO_export_univId"), 
                "Save Background Gene ID's to file", 
                "Save Background Gene ID's to file...", 
                filetype = list(txt = "txt")),
        ),
        column(9,
            plotlyOutput(ns("plot_GO"), height = "800px")
        )
    )
}