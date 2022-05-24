ui_sidebar <- function() {
    dashboardSidebar(
        sidebarMenu(id = "navSelection",
            menuItem("About", tabName = "navTitle"),
            # menuItem("System", tabName = "navSystem"),
            menuItem("Reference", tabName = "navRef_New"),
            menuItem("Experiment", 
                tabName = "navExpr"),
            menuItem("Analysis",
                menuSubItem("Load Experiment", tabName = "navExprLoad"),
                menuSubItem("Experiment QC", tabName = "navQC"),
                menuSubItem("Filters", tabName = "navFilter"),
                menuSubItem("Differential Expression Analysis", 
                    tabName = "navAnalyse")
            ),
            menuItem("Display",
                menuSubItem("Volcano Plot", tabName = "navVolcano"),
                menuSubItem("Scatter Plot", tabName = "navDiag"),
                menuSubItem("Heatmap", tabName = "navHeatmap"),
                menuSubItem("Coverage Plot", tabName = "navCoverage")
            )
        )
    )
}

ui_tab_title <- function() {
    max_threads <- parallel::detectCores()
    thread_choices <- c(1, 2, 4, 6, 8, 12, 16, 24, 32)
    thread_choices <- thread_choices[thread_choices <= max_threads]
    out_choices <- c(thread_choices, "custom")
    tabItem(tabName = "navTitle",
        fluidRow(
            column(4,
            tags$div(title = paste(
                "These options set the number of threads to run",
                "computationally-intensive operations",
                "such as processBAM, collateData, and DESeq2.",
                "Also, set memory mode here."),
                h3("Thread and memory options"),
                shinyWidgets::pickerInput(
                   inputId = "thread_number",
                   label = "Number of threads to use", 
                    choices = out_choices
                ),
                conditionalPanel(
                    condition = "['custom'].indexOf(input.thread_number) >= 0",
                    numericInput("cores_numeric", "# Threads", min = 1, 
                        max = parallel::detectCores(), value = 1)
                ),
                radioGroupButtons(
                    inputId = "memory_option",
                    label = "Memory Usage",
                    choices = c( 
                        "Low",
                        "High"
                    ),
                    justified = TRUE,
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon")
                    )
                ), br(),
                h4("Estimated memory usage for reference generation"),
                textOutput("txt_mem_buildRef"), br(),
                h4("Estimated memory usage for BAM processing"),
                textOutput("txt_mem_processBAM"), br(),
                h4("Estimated memory usage for data collation"),
                textOutput("txt_mem_collateData"), br()
            )),
            column(8,
                img(src=paste0(
                    "https://pbs.twimg.com/",
                    "profile_images/",
                    "1310789966293655553/",
                    "7HawCItY_400x400.jpg"
                ))
            )
        )
    )
}

ui_tab_ref_new <- function() {
    tabItem(tabName = "navRef_New",
        ui_ref_new("new_ref")
    )
}

ui_tab_expr <- function() {
    tabItem(tabName = "navExpr",
        ui_expr("build_expr")
    )
}

ui_tab_expr_load <- function() {
    tabItem(tabName = "navExprLoad",
        ui_expr_limited("load_expr")
    )
}

ui_tab_qc <- function() {
    tabItem(tabName = "navQC",
        ui_qc("qc")
    )
}

ui_tab_filter <- function() {
    tabItem(tabName = "navFilter",
        # Current Experiment
        ui_filters("filters")
    )
}

ui_tab_analyse <- function() {
    tabItem(tabName = "navAnalyse",
        ui_DE("DE")
    )
}

ui_tab_diag <- function() {
    tabItem(tabName = "navDiag",
        ui_vis_diag("diag")
    )
}

ui_tab_volcano <- function() {
    tabItem(tabName = "navVolcano",
        ui_vis_volcano("volcano")
    )
}

ui_tab_heatmap <- function() {
    tabItem(tabName = "navHeatmap",
        ui_vis_heatmap("heatmap")
    )
}

ui_tab_coverage <- function() {
    tabItem(tabName = "navCoverage",
        ui_cov("cov")
    )
}

################################################################################

ui_toggle_wellPanel <- function(
        inputId, title, color = "danger", icon = icon("bars"), ...
) {
    tagList(
        actionBttn(
            inputId = inputId,
            label = title,
            style = "gradient", 
            color = color,
            icon = icon
        ),
        br(),
        conditionalPanel(condition = 
            paste0("input.", inputId, " % 2 != 0"),
            wellPanel(...)
        )
    )
}

ui_toggle_wellPanel_modular <- function(
        inputId, id, title, color = "danger", icon = icon("bars"), ...
) {
    ns <- NS(id)
    tagList(
        actionBttn(
            inputId = ns(inputId),
            label = title,
            style = "gradient", 
            color = color,
            icon = icon
        ),
        br(),
        conditionalPanel(
            ns = ns,
            condition = paste0("input.", inputId, " % 2 != 0"),
            wellPanel(...)
        )
    )
}



