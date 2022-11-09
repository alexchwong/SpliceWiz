ui_DE <- function(id) {
    ns <- NS(id)
    DE_opts <- c()
    if(.check_package_installed("limma", "3.44.0", "silent")) {
        DE_opts <- "limma"
    }
    if(.check_package_installed("DESeq2", "1.30.0", "silent")) {
        DE_opts <- c(DE_opts, "DESeq2")
    }
    if(.check_package_installed("DoubleExpSeq", "1.1", "silent")) {
        DE_opts <- c(DE_opts, "DoubleExpSeq")
    }
    if(.check_package_installed("satuRn", "1.4.2", "silent")) {
        DE_opts <- c(DE_opts, "satuRn")
    }
    if(is.null(DE_opts)) {
        DE_opts <- ("(none) - please install limma, DESeq2, DoubleExpSeq, or satuRn")
    }
    fluidRow(
        .ui_notice(),
        column(4,
            textOutput(ns("warning_DE")),
            selectInput(ns('method_DE'), 'Method', 
                DE_opts),
            selectInput(ns('variable_DE'), 'Variable', 
                c("(none)")),
            selectInput(ns('nom_DE'), 'Nominator', 
                c("(none)")),
            selectInput(ns('denom_DE'), 'Denominator', 
                c("(none)")),
            selectInput(ns('batch1_DE'), 'Batch Factor 1', 
                c("(none)")),
            selectInput(ns('batch2_DE'), 'Batch Factor 2', 
                c("(none)")),
            # conditionalPanel(ns = ns,
                # condition = "['(time series)'].indexOf(input.nom_DE) >= 0",
                # numericInput(ns("degrees_DE"), 
                    # label = "Splines degrees of freedom (limma only)", 
                    # value = 1, min = 1, max = 5)
            # ),
            shinyWidgets::radioGroupButtons(ns("de_IRmode"), 
                label = "Intron Retention analysis (mode)",
                    justified = FALSE,
                choiceNames = c(
                    "All Introns (IR-Ratio)", 
                    "Annotated IR events (IR-Ratio)", 
                    "Annotated IR events (binary PSI)"
                ), 
                choiceValues = c("all", "annotated", "annotated_binary"),
                checkIcon = list(yes = icon("ok", lib = "glyphicon"))
            ),
            shinyWidgets::switchInput(ns("adjP_DE"), 
                label = "Sort by Adjusted P Values", 
                value = TRUE, labelWidth = "120px"),
            actionButton(ns("perform_DE"), "Perform DE"),
            # actionButton(ns("perform_DE_spoof"), "Perform DE (spoof)"),
            shinyFilesButton(ns("load_DE"), label = "Load DE", 
                    title = "Load DE from RDS", multiple = FALSE),
            shinySaveButton(ns("save_DE"), "Save DE", "Save DE as...", 
                filetype = list(RDS = "Rds")),
        ),
        column(8,
            actionButton(ns("clear_selected_DE"), "Clear Selected Events"),
            div(style = 'overflow-x: scroll',  
                DT::dataTableOutput(ns('DT_DE'))
            )
        )
    )
}