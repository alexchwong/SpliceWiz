ui_GO <- function(id) {
    ns <- NS(id)
    fluidRow(
        .ui_notice(),
        column(3,
            visFilter_UI(ns("GOfilters")),
            textOutput(ns("warning_GO")),
            selectInput(ns('category_GO'), 'GO Category', 
                c("Biological Pathway", "Molecular Function", 
                    "Cellular Compartment")),
            shinyWidgets::radioGroupButtons(ns("direction_GO"), 
                label = "Regulation direction", justified = FALSE,
                choices = c("Up+Down", "Up", "Down"), 
                checkIcon = list(yes = icon("ok", lib = "glyphicon"))
            ),
            selectInput(ns('universe_GO'), 'Background Genes based on', 
                c("All ASE events", "Selected ASE Modality", "All Genes")),
            conditionalPanel(ns = ns,
                condition = paste0(
                    "['Selected ASE Modality'].",
                    "indexOf(input.universe_GO) == 0"
                ),
                selectInput(ns("GO_EventType"), 
                    "Filter Events by ASE Modality", 
                    width = '100%', multiple = TRUE,
                    choices = c("IR", "MXE", "SE", "AFE", "ALE", 
                        "A5SS", "A3SS")), br(),
            ),
            vis_ggplot_UI(ns("GOplotSave")),
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