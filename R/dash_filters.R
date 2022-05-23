ui_filters <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(4,
                wellPanel(style = "overflow-y:scroll; max-height: 800px",
                    filterModule_UI(ns("filter1"), "Filter #1"),
                    filterModule_UI(ns("filter2"), "Filter #2"),
                    filterModule_UI(ns("filter3"), "Filter #3"),
                    filterModule_UI(ns("filter4"), "Filter #4"),
                    filterModule_UI(ns("filter5"), "Filter #5"),
                    filterModule_UI(ns("filter6"), "Filter #6")
                )
            ),
            column(4,
                wellPanel(style = "overflow-y:scroll; max-height: 800px",
                    filterModule_UI(ns("filter7"), "Filter #7"),
                    filterModule_UI(ns("filter8"), "Filter #8"),
                    filterModule_UI(ns("filter9"), "Filter #9"),
                    filterModule_UI(ns("filter10"), "Filter #10"),
                    filterModule_UI(ns("filter11"), "Filter #11"),
                    filterModule_UI(ns("filter12"), "Filter #12")
                )
            ),
            column(4,
                textOutput(ns("current_expr_Filters")), br(),
                textOutput(ns("current_ref_Filters")), br(),
                actionButton(ns("loadDefault_Filters"), "Load Default Filters"),
                br(), br(),
                actionButton(ns("refresh_filters_Filters"), "Apply Filters"),
                br(), br(),
                shinySaveButton(ns("saveAnalysis_Filters"), "Save Filters", 
                    "Save Filters as...", filetype = list(RDS = "Rds")),
                shinyFilesButton(ns("loadAnalysis_Filters"), 
                    label = "Load Filters", 
                    title = "Load Filters from Rds", multiple = FALSE),
                plotlyOutput(ns("plot_filtered_Events")),
                selectInput(ns('graphscale_Filters'), 'Y-axis Scale', 
                    width = '100%', choices = c("linear", "log10")), 
            )
        )
    )
}

server_filters <- function(
        id, refresh_tab, volumes, get_se, get_filters_from_DE
) {
    moduleServer(id, function(input, output, session) {
        settings_filter <- setreactive_filtered_SE()

        # Tab refresh
        observeEvent(refresh_tab(), {
            if(is_valid(get_se())) {
                output$current_expr_Filters <- 
                    renderText("NxtSE loaded")
                processFilters()
            } else {
                output$current_expr_Filters <- 
                    renderText("Please load NxtSE first")
            }
        })

        # Reactives to individual filter modules
        getFilterData <- function(i) {
            if(
                is_valid(settings_filter$filters) &&
                length(settings_filter$filters) >= i
            
            ) {
                return(settings_filter$filters[[i]])
            } else {
                return(list())
            }
        }
        r_filter1 <- reactive({getFilterData(1)})
        r_filter2 <- reactive({getFilterData(2)})
        r_filter3 <- reactive({getFilterData(3)})
        r_filter4 <- reactive({getFilterData(4)})
        r_filter5 <- reactive({getFilterData(5)})
        r_filter6 <- reactive({getFilterData(6)})
        r_filter7 <- reactive({getFilterData(7)})
        r_filter8 <- reactive({getFilterData(8)})
        r_filter9 <- reactive({getFilterData(9)})
        r_filter10 <- reactive({getFilterData(10)})
        r_filter11 <- reactive({getFilterData(11)})
        r_filter12 <- reactive({getFilterData(12)})
        
        conditionList <- reactive({
            req(get_se())
            if(is(get_se(), "NxtSE")) {
                colnames(colData(get_se()))
            } else {
                c("")
            }
        })
        filter1 <- filterModule_server("filter1", r_filter1, conditionList)
        filter2 <- filterModule_server("filter2", r_filter2, conditionList)
        filter3 <- filterModule_server("filter3", r_filter3, conditionList)
        filter4 <- filterModule_server("filter4", r_filter4, conditionList)
        filter5 <- filterModule_server("filter5", r_filter5, conditionList)
        filter6 <- filterModule_server("filter6", r_filter6, conditionList)
        filter7 <- filterModule_server("filter7", r_filter7, conditionList)
        filter8 <- filterModule_server("filter8", r_filter8, conditionList)
        filter9 <- filterModule_server("filter9", r_filter9, conditionList)
        filter10 <- filterModule_server("filter10", r_filter10, conditionList)
        filter11 <- filterModule_server("filter11", r_filter11, conditionList)
        filter12 <- filterModule_server("filter12", r_filter12, conditionList)


        # Function to process NxtSE with set filters
        processFilters <- function() {
            message("Refreshing filters")
            if(is(get_se(), "NxtSE")) {
                filterSummary <- rep(TRUE, nrow(get_se()))
                if(is_valid(settings_filter$filters)) {
                    filters_to_run <- c()
                    for(i in seq_len(12)) {
                        if(
                            length(settings_filter$filters) >= i &&
                            is_valid(settings_filter$filters[[i]]@filterType)  
                        ) {
                            filters_to_run <- c(filters_to_run, i)
                        }
                    }
                    withProgress(message = 'Running NxtSE Filters', value = 0, {
                        for(i in filters_to_run) {
                            filterSummary <- filterSummary & runFilter(
                                get_se(),
                                settings_filter$filters[[i]]
                            )
                            incProgress(1/length(filters_to_run))
                        }
                    })
                    .filters_sweetalert_finish(session, 
                        sum(filterSummary == TRUE))
                }
                settings_filter$filterSummary <- filterSummary
                message("Filtered ", sum(filterSummary == TRUE), " ASE events")
            } else if(is_valid(settings_filter$filters)) {
                for(i in seq_len(12)) {
                    if(
                        length(settings_filter$filters) >= i &&
                        is_valid(settings_filter$filters[[i]]@filterType)  
                    ) {
                        print(settings_filter$filters[[i]])
                    }
                }
            }
        }
            
        # Function to run filters using above function
        observeEvent({list(
            input$refresh_filters_Filters
        )}, {
            req(input$refresh_filters_Filters)
            settings_filter$filters[[1]] <- filter1()
            settings_filter$filters[[2]] <- filter2()
            settings_filter$filters[[3]] <- filter3()
            settings_filter$filters[[4]] <- filter4()
            settings_filter$filters[[5]] <- filter5()
            settings_filter$filters[[6]] <- filter6()
            settings_filter$filters[[7]] <- filter7()
            settings_filter$filters[[8]] <- filter8()
            settings_filter$filters[[9]] <- filter9()
            settings_filter$filters[[10]] <- filter10()
            settings_filter$filters[[11]] <- filter11()
            settings_filter$filters[[12]] <- filter12()
            processFilters()
        })
        
        # Updates column chart based on changing applied filter, or log scale
        observeEvent({list(
            settings_filter$filterSummary,
            input$graphscale_Filters
        )}, {
            req(get_se())
            req(settings_filter$filterSummary)
            DT <- data.table(
                EventType = rowData(get_se())$EventType,
                keep = settings_filter$filterSummary
            )
            output$plot_filtered_Events <- renderPlotly({
                print(
                    Filters_Plot_Summary(DT, input$graphscale_Filters)
                )
            })
        })
        
        # Saves current list of filters to Rds file
        observe({
            shinyFileSave(input, "saveAnalysis_Filters", 
                roots = volumes(), session = session)        
        })
        observeEvent(input$saveAnalysis_Filters, {
            req(settings_filter$filters)
            selectedfile <- parseSavePath(volumes(), input$saveAnalysis_Filters)
            req(selectedfile$datapath)

            settings_filter$filters[[1]] <- filter1()
            settings_filter$filters[[2]] <- filter2()
            settings_filter$filters[[3]] <- filter3()
            settings_filter$filters[[4]] <- filter4()
            settings_filter$filters[[5]] <- filter5()
            settings_filter$filters[[6]] <- filter6()
            settings_filter$filters[[7]] <- filter7()
            settings_filter$filters[[8]] <- filter8()
            settings_filter$filters[[9]] <- filter9()
            settings_filter$filters[[10]] <- filter10()
            settings_filter$filters[[11]] <- filter11()
            settings_filter$filters[[12]] <- filter12()
            
            final <- settings_filter$filters
            saveRDS(final, selectedfile$datapath)
        })

        # Load filter list from Rds file
        observe({
            shinyFileChoose(input, "loadAnalysis_Filters", 
                roots = volumes(), session = session,
                filetypes = c("Rds"))        
        })
        observeEvent(input$loadAnalysis_Filters, {
            selectedfile <- parseFilePaths(volumes(), 
                input$loadAnalysis_Filters)
            req(selectedfile$datapath)
            settings_filter$filters <- readRDS(selectedfile$datapath)
        })
        
        # Import filters from loading DE object
        observeEvent(get_filters_from_DE(), {
            req(get_filters_from_DE())
            settings_filter$filters = get_filters_from_DE()
        })
        observeEvent(input$loadDefault_Filters, {
            settings_filter$filters = getDefaultFilters()
        })
        return(settings_filter)
    })
}

# Helper function to plot the number of filtered ASE
Filters_Plot_Summary <- function(DT, scale) {
    if(scale == "log10") {
        DT[, c("Included", "Excluded") := .(
            log10(sum(get("keep") == TRUE)),
            log10(sum(!is.na(get("keep")))) - 
                log10(sum(get("keep") == TRUE))
        ), by = "EventType"]
    } else {
        DT[, c("Included", "Excluded") := .(
            sum(get("keep") == TRUE),
            sum(get("keep") != TRUE)
        ), by = "EventType"] 
    }
    DT <- unique(DT, by = "EventType")
    incl <- as.data.frame(DT[, c("EventType", "Included")])
    incl$filtered <- "Included"
    colnames(incl)[colnames(incl) == "Included"] <- "Events"
    excl <- as.data.frame(DT[, c("EventType", "Excluded")])
    excl$filtered <- "Excluded"
    colnames(excl)[colnames(excl) == "Excluded"] <- "Events"    

    # ggplot summary as bar plot
    p <- ggplot(rbind(incl, excl), 
        aes(x = get("EventType"), y = get("Events"), fill = get("filtered"))
    ) + geom_bar(position="stack", stat="identity")
    if(scale == "log10") {
        p <- p + labs(x = "Event Type", y = "log10 Events")
    } else {
        p <- p + labs(x = "Event Type", y = "Events")        
    }
    p <- p + labs(fill = "Filtered")
    return(p)
}

.filters_sweetalert_finish <- function(session, num_events) {
    sendSweetAlert(
        session = session,
        title = paste(
            "NxtSE filters processed (",
            num_events, " ASEs retained)"
        ),
        type = "success"
    )
}