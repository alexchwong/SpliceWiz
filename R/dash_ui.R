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
                menuSubItem("Gene Ontology Analysis", tabName = "navGO"),
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
            column(8,
                h2("Welcome to SpliceWiz!"),
                # tags$div(title = paste(
                    # "These options set the number of threads to run",
                    # "computationally-intensive operations",
                    # "such as processBAM, collateData, and DESeq2.",
                    # "Also, set memory mode here."),
                    # h3("Thread and memory options"),
                    # shinyWidgets::pickerInput(
                       # inputId = "thread_number",
                       # label = "Number of threads to use", 
                        # choices = out_choices
                    # ),
                    # conditionalPanel(
                        # condition = "['custom'].indexOf(input.thread_number) >= 0",
                        # numericInput("cores_numeric", "# Threads", min = 1, 
                            # max = parallel::detectCores(), value = 1)
                    # ),
                    # radioGroupButtons(
                        # inputId = "memory_option",
                        # label = "Memory Usage",
                        # choices = c( 
                            # "Low",
                            # "High"
                        # ),
                        # justified = TRUE,
                        # checkIcon = list(
                            # yes = icon("ok", lib = "glyphicon")
                        # )
                    # ), br(),
                    # h4("Estimated memory usage for reference generation"),
                    # textOutput("txt_mem_buildRef"), br(),
                    # h4("Estimated memory usage for BAM processing"),
                    # textOutput("txt_mem_processBAM"), br(),
                    # h4("Estimated memory usage for data collation"),
                    # textOutput("txt_mem_collateData"), br()
                # ),
                img(src=paste0(
                    "https://www.biorxiv.org/content/biorxiv/", 
                    "early/2022/07/06/2022.07.05.498887/", 
                    "F1.large.jpg"                
                ), height = 800)                
            ),
            column(4,
                # img(src=paste0(
                    # "https://pbs.twimg.com/",
                    # "profile_images/",
                    # "1310789966293655553/",
                    # "7HawCItY_400x400.jpg"
                # ))
                img(src="localImages/labLogo.jpg")
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

ui_tab_GO <- function() {
    tabItem(tabName = "navGO",
        ui_GO("GO")
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
        # ui_cov("cov")
        ui_cov_new("cov")
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

shinyDirAttnButton <- function(
    id, label, title,
    icon = NULL,
    style = "gradient",
    color = "default",
    size = "md",
    block = FALSE,
    no_outline = TRUE,
  ...
) {
    value <- restoreInput(id = id, default = NULL)
    style <- match.arg(
        arg = style,
        choices = c(
            "simple", "bordered", "minimal", "stretch", "jelly",
            "gradient", "fill", "material-circle", "material-flat",
            "pill", "float", "unite"
        )
    )
    color <- match.arg(
        arg = color,
        choices = c(
            "default", "primary", "warning", "danger", "success", "royal"
        )
    )
    size <- match.arg(arg = size, choices = c("xs", "sm", "md", "lg"))
  
    tagList(
        singleton(
            tags$head(
                tags$script(src = "sF/shinyFiles.js"),
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "sF/styles.css"
                ),
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "sF/fileIcons.css"
                )
            )
        ),
        sWidgets_attachDep(
            tags$button(
                id = id,
                type = "button",
                class = "shinyDirectories action-button bttn",
                class = paste0("bttn-", style),
                class = paste0("bttn-", size),
                class = paste0("bttn-", color),
                class = if (block) "bttn-block",
                class = if (no_outline) "bttn-no-outline",
                style = style,
                "data-title" = title,
                "data-val" = value,
                icon, label, ...
            ), "bttn"
        )
    )
}

shinyFilesAttnButton <- function(
    id, label, title, 
    multiple = FALSE, 
    buttonType="default", 
    icon = NULL,
    style = "gradient",
    color = "default",
    size = "md",
    block = FALSE,
    no_outline = TRUE,
    viewtype = c("detail", "list", "icon"), 
    ...
) {
    value <- restoreInput(id = id, default = NULL)
    viewtype <- match.arg(viewtype)
    if(!is_valid(viewtype)) viewtype <- "detail"
    style <- match.arg(
        arg = style,
        choices = c(
            "simple", "bordered", "minimal", "stretch", "jelly",
            "gradient", "fill", "material-circle", "material-flat",
            "pill", "float", "unite"
        )
    )
    color <- match.arg(
        arg = color,
        choices = c(
            "default", "primary", "warning", "danger", "success", "royal"
        )
    )
    size <- match.arg(arg = size, choices = c("xs", "sm", "md", "lg"))
    
    tagList(
        singleton(
            tags$head(
                tags$script(src = "sF/shinyFiles.js"),
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "sF/styles.css"
                ),
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "sF/fileIcons.css"
                )
            )
        ),
        sWidgets_attachDep(
            tags$button(
                id = id,
                type = "button",
                class = "shinyFiles action-button bttn",
                class = paste0("bttn-", style),
                class = paste0("bttn-", size),
                class = paste0("bttn-", color),
                class = if (block) "bttn-block",
                class = if (no_outline) "bttn-no-outline",
                style = style,
                "data-title" = title,
                "data-selecttype" = ifelse(multiple, "multiple", "single"),
                "data-val" = value,
                "data-view" = paste0("sF-btn-", viewtype),
                icon, label, ...
            ), "bttn"
        )
    )
}

shinySaveAttnButton <- function(
    id, label, title, 
    filename = "", filetype,
    buttonType="default", 
    icon = NULL,
    style = "gradient",
    color = "default",
    size = "md",
    block = FALSE,
    no_outline = TRUE,
    viewtype = c("detail", "list", "icon"), 
    ...
) {
    value <- restoreInput(id = id, default = NULL)
    viewtype <- match.arg(viewtype)
    if (missing(filetype)) filetype <- NA
    filetype <- sFiles_formatFiletype(filetype)
    if(!is_valid(viewtype)) viewtype <- "detail"
    style <- match.arg(
        arg = style,
        choices = c(
            "simple", "bordered", "minimal", "stretch", "jelly",
            "gradient", "fill", "material-circle", "material-flat",
            "pill", "float", "unite"
        )
    )
    color <- match.arg(
        arg = color,
        choices = c(
            "default", "primary", "warning", "danger", "success", "royal"
        )
    )
    size <- match.arg(arg = size, choices = c("xs", "sm", "md", "lg"))
    
    tagList(
        singleton(
            tags$head(
                tags$script(src = "sF/shinyFiles.js"),
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "sF/styles.css"
                ),
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "sF/fileIcons.css"
                )
            )
        ),
        sWidgets_attachDep(
            tags$button(
                id = id,
                type = "button",
                class = "shinySave action-button bttn",
                class = paste0("bttn-", style),
                class = paste0("bttn-", size),
                class = paste0("bttn-", color),
                class = if (block) "bttn-block",
                class = if (no_outline) "bttn-no-outline",
                style = style,
                "data-title" = title,
                "data-filetype" = filetype,
                "data-filename" = filename,
                "data-val" = value,
                "data-view" = paste0("sF-btn-", viewtype),
                icon, label, ...
            ), "bttn"
        )
    )
}
