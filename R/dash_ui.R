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
                menuSubItem("Scatter Plot", tabName = "navDiag"),
                menuSubItem("Volcano Plot", tabName = "navVolcano"),
                menuSubItem("Heatmap", tabName = "navHeatmap"),
                menuSubItem("Coverage Plot", tabName = "navCoverage")
            )
        )
    )
}

ui_tab_title <- function() {
    tabItem(tabName = "navTitle",
        box(
            tags$div(title = paste(
                "Number of threads to run",
                "computationally-intensive operations",
                "such as processBAM, collateData, and DESeq2"),
                radioGroupButtons(
                    inputId = "thread_option",
                    label = "Mode",
                    choices = c( 
                        "Multi-Thread (High)",
                        "Multi-Thread (Low)", 
                        "Single-Thread", "Custom"
                    ),
                    justified = TRUE,
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon")
                    )
                ),
                conditionalPanel(
                    condition = "['Custom'].indexOf(input.thread_option) >= 0",
                    numericInput("cores_numeric", "# Threads", min = 1, 
                        max = parallel::detectCores(), value = 1)
                )                
            )
        ), br(),
        img(src=paste0(
            "https://pbs.twimg.com/",
            "profile_images/",
            "1310789966293655553/",
            "7HawCItY_400x400.jpg"
        ))
    )
}

ui_tab_system <- function() {
    tabItem(tabName = "navSystem",
        box(
            tags$div(title = paste(
                "Number of threads to run",
                "computationally-intensive operations",
                "such as processBAM, collateData, and DESeq2"),
                radioGroupButtons(
                    inputId = "thread_option",
                    label = "Mode",
                    choices = c( 
                        "Multi-Thread (High)",
                        "Multi-Thread (Low)", 
                        "Single-Thread", "Custom"
                    ),
                    justified = TRUE,
                    checkIcon = list(
                        yes = icon("ok", lib = "glyphicon")
                    )
                ),
                conditionalPanel(
                    condition = "['Custom'].indexOf(input.thread_option) >= 0",
                    numericInput("cores_numeric", "# Threads", min = 1, 
                        max = parallel::detectCores(), value = 1)
                )                
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
    ns = NS(id)
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



