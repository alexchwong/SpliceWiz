# Order of DE filtration
# get_de()
# - get_de()[rows_all]
# - get_de()[rows_all][seq_len(Top N rows)]

# New filtration approach
# get_de()
# - get_de()[rows_all]
#   - {option to filter by padj, pvalue, or top n rows}

vis_ggplot_server <- function(id, get_ggplot, volumes) {
    moduleServer(id, function(input, output, session) {
        final <- reactiveVal(value = NULL)
        observe({
            shinyFileSave(input, "saveplot", roots = volumes(), 
                session = session, filetypes = c("pdf"))
        })
        observeEvent(input$saveplot, {
            p <- get_ggplot()
            req(p)

            selectedfile <- parseSavePath(volumes(), input$saveplot)
            req(selectedfile$datapath)
            
            ggsave(
                selectedfile$datapath, 
                plot = p,
                width = isolate(input$wt), 
                height = isolate(input$ht), 
                units = isolate(input$units)
            )
        })
        return(final)
    })
}

visFilter_server <- function(id, get_de, rows_all, rows_selected) {
    moduleServer(id, function(input, output, session) {
        final <- reactiveVal(
            value = data.table()
        )
        
        # observeEvent(input$vF_reset, {
            # updateSelectInput(inputId = "vF_filterType", 
                # selected = "Adjusted P value")
            # updateSliderInput(inputId = "vF_pvalT", value = 0.05)
            # updateSliderInput(inputId = "vF_topN", value = 500)
            # updateSelectInput(inputId = "vF_EventType", selected = NULL)
        # })
        
        input_df_r <- reactive({
            req(get_de())
            res_bkgd <- as.data.table(.get_unified_volcano_data(get_de()))
            res_bkgd$All <- seq_len(nrow(res_bkgd)) %in% rows_all()
            res_bkgd$Selected <- seq_len(nrow(res_bkgd)) %in% rows_selected()
            return(res_bkgd)
        })
        
        useDE_r <- reactive({
            req(input_df_r())
            tmpres <- input_df_r()
            # message("Original #events: ", nrow(tmpres))
            if(is_valid(input$vF_EventType)) {
                tmpres <- tmpres[get("EventType") %in% input$vF_EventType]
                # message("after Event filter: ", nrow(tmpres))
            }

            req(input$vF_filterType)
            if(input$vF_filterType == "Adjusted P value") {
                tmpres <- tmpres[get("FDR") <= input$vF_pvalT]
            } else if(input$vF_filterType == "Nominal P value") {
                tmpres <- tmpres[get("pvalue") <= input$vF_pvalT]
            } else if(input$vF_filterType == "Top events by p-value") {
                if(input$vF_topN < nrow(tmpres)) {
                    tmpres <- tmpres[seq_len(input$vF_topN)]
                }
            } else if(input$vF_filterType == "Highlighted events") {
                tmpres <- tmpres[get("Selected") == TRUE]
            }
            # message("after final filter: ", nrow(tmpres))
            return(tmpres)
        })
        observe(final(useDE_r()))
        return(final)
    })
}


server_vis_diag <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Diag <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
            randomVar <- runif(1)
            settings_Diag$trigger <- randomVar
            # print(randomVar)
        })
        
        observeEvent(settings_Diag$trigger, {
            req(settings_Diag$trigger)
            
            output$warning_diag <- renderText({
                validate(need(is(get_se(), "NxtSE"), 
                "Please load NxtSE object"))
                
                "NxtSE Loaded"
            })

            # Update annotation column names in selection
            req(is(get_se(), "NxtSE"))
            colData <- colData(get_se())
            if(
                    is_valid(input$variable_diag) && 
                    input$variable_diag %in% colnames(colData)
            ) {
                selectedOption <- isolate(input$variable_diag)
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = selectedOption
                )
            } else {
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = "(none)"
                )
            }
        })
        
        useDE_r <- visFilter_server("scatter", get_de, rows_all, rows_selected)

        # Reactive to generate mean PSIs
        # - this is the main bottleneck
        observe({
            req(settings_Diag$trigger)
            req(is(get_se(), "NxtSE"))
            req(get_de())
            
            req(is_valid(input$variable_diag))
            req(is_valid(input$denom_diag))
            req(is_valid(input$nom_diag))

            tmpres <- as.data.table(
                .get_unified_volcano_data(get_de()[rows_all(),]))

            req(tmpres$EventName)
            req(all(tmpres$EventName %in% rowData(get_se())$EventName))

            withProgress(message = 'Calculating mean PSIs...', value = 0, {
                df.diag <- makeMeanPSI(
                    get_se(), tmpres$EventName, input$variable_diag, 
                    list(input$nom_diag, input$denom_diag)
                )
                colnames(df.diag)[seq(2,3)] <- c("nom", "denom")
                settings_Diag$meanPSI <- df.diag
            })
        })

        # Update local rows_selected with that of global
        observeEvent(rows_selected(), {
            settings_Diag$selected <- rows_selected()
        }, ignoreNULL = FALSE)
    
        # Main scatter plot generator
        observe({
            req(settings_Diag$trigger)
            req(is(get_se(), "NxtSE"))
            req(get_de())

            req(settings_Diag$meanPSI)
            
            res <- useDE_r()
            req(res)
            req(nrow(res) > 0)
            req(all(res$EventName %in% rownames(get_se())))

            df.diag <- res[get("EventName") %in% 
                settings_Diag$meanPSI$EventName]
            req(nrow(df.diag) > 0)
            
            if(input$NMD_diag == TRUE) {
                df.diag <- df.diag[df.diag$NMD_direction != 0, ]
            }
            req(nrow(df.diag) > 0)
            
            df.diag$nom <- settings_Diag$meanPSI$nom[match(
                df.diag$EventName, settings_Diag$meanPSI$EventName
            )]
            df.diag$denom <- settings_Diag$meanPSI$denom[match(
                df.diag$EventName, settings_Diag$meanPSI$EventName
            )]
                        
            # Annotate which rows are selected; NMD direction
            selected <- settings_Diag$selected      
            if(is_valid(selected)) {
                df.diag$selected <- 
                    (df.diag$EventName %in% get_de()$EventName[selected])
            } else {
                df.diag$selected <- FALSE
            }
            df.diag$NMD_direction <- .getNMDcode(get_de()$flags[
                match(df.diag$EventName, get_de()$EventName)])
            
            if(input$NMD_diag == TRUE) {
                df.diag$nom_NMD     <- ifelse(df.diag$NMD_direction == 1, 
                                        df.diag$nom, df.diag$denom)
                df.diag$denom_NMD   <- ifelse(df.diag$NMD_direction == -1, 
                                        df.diag$nom, df.diag$denom)
            }

            settings_Diag$plot_ini <- TRUE            
            if(input$NMD_diag == TRUE) {
                p <- ggplot(df.diag, 
                        aes(
                            x = get("nom_NMD"), y = get("denom_NMD"), 
                            key = get("EventName"), text = get("EventName"), 
                            colour = get("selected")
                        )
                    ) + geom_point() + 
                    scale_color_manual(values = c("black", "red")) +
                    labs(
                        x = paste(input$nom_diag, "NMD substrate"),
                        y = paste(input$denom_diag, "NMD substrate")
                    )
            } else {
                p <- ggplot(df.diag, 
                        aes(
                            x = get("nom"), y = get("denom"), 
                            key = get("EventName"), text = get("EventName"), 
                            colour = get("selected")
                        )
                    ) + geom_point() + 
                    scale_color_manual(values = c("black", "red")) +
                    labs(
                        x = paste(input$nom_diag),
                        y = paste(input$denom_diag)
                    )         
            }
            
            # Annotate colors, etc
            p <- p + labs(color = "Selected")
            settings_Diag$ggplot <- p
            
            withProgress(message = 'Rendering plot...', value = 0, {
                # Record ggplot / plotly objects into settings_Diag
                py <- ggplotly(
                    p, tooltip = "text",
                    source = "plotly_diagonal",
                    type = "scatter_gl"
                ) %>% toWebGL() %>% layout(
                    dragmode = "lasso",
                    yaxis = list(scaleanchor="x", scaleratio=1)
                )
                # Add hoveron entry
                py$x$data <- lapply(py$x$data, function(x) {
                    x$hoveron <- NULL
                    x
                })
                settings_Diag$final_plot <- py
                
                if(packageVersion("plotly") >= "4.9.0") {
                    plotly::event_register(
                        settings_Diag$final_plot, "plotly_click")
                    plotly::event_register(
                        settings_Diag$final_plot, "plotly_selected")
                }
            })
        })
    
        output$plot_diag <- renderPlotly({
            validate(need(is(get_se(), "NxtSE"), "Load Experiment first"))
            validate(need(get_de(), "Perform DE Analysis first"))

            validate(need(is_valid(input$variable_diag), 
                "Select conditions and contrasts"))
            validate(need(is_valid(input$denom_diag), 
                "Select conditions and contrasts"))
            validate(need(is_valid(input$nom_diag), 
                "Select conditions and contrasts"))
        
            # Filter DE by EventType; fetch diag object
            res <- useDE_r()
            validate(need(nrow(res) > 0, 
                "No events found. Consider relaxing some filters"
            ))
            validate(need(all(res$EventName %in% rownames(get_se())),
                "Filtered NxtSE does not match DE results"
            ))
            
            df.diag <- res[get("EventName") %in% 
                settings_Diag$meanPSI$EventName]
            validate(need(nrow(df.diag) > 0, 
                "No events found. Consider relaxing some filters"
            ))

            if(input$NMD_diag == TRUE) {
                df.diag <- df.diag[df.diag$NMD_direction != 0, ]
            }
            validate(need(nrow(df.diag) > 0, 
                paste(
                    "No events found.", 
                    "Consider relaxing some filters or disable NMD mode"
                )
            ))
            
            req(settings_Diag$final_plot)
            settings_Diag$final_plot
        })
        
        # Reactive click
        settings_Diag$plotly_click <- reactive({
            plot_exist <- settings_Diag$plot_ini
            if(plot_exist) 
                event_data("plotly_click", source = "plotly_diagonal")
        })
        observeEvent(settings_Diag$plotly_click(), {
            req(settings_Diag$plotly_click())
            click <- settings_Diag$plotly_click()
            click.id <- which(get_de()$EventName == click$key)
            req(click.id)

            selected <- settings_Diag$selected

            if(click.id %in% selected && input$reverse_select) {
                selected <- selected[-which(selected == click.id)]
            } else if(!input$reverse_select) {
                selected <- c(selected, click.id)
            }
            settings_Diag$selected <- selected
        })

        # Reactive brush
        settings_Diag$plotly_brush <- reactive({
            plot_exist <- settings_Diag$plot_ini
            if(plot_exist)
                event_data("plotly_selected", source = "plotly_diagonal")
        })
        observeEvent(settings_Diag$plotly_brush(), {
            req(settings_Diag$plotly_brush())
            brush <- settings_Diag$plotly_brush()
            brush.id <- which(get_de()$EventName %in% brush$key)
            req(brush.id)

            selected <- settings_Diag$selected
            if(!input$reverse_select) {
                selected <- union(selected, brush.id)           
            } else {
                selected <- setdiff(selected, brush.id)
            }
            settings_Diag$selected <- selected
        })
    
        # Update nominator / denominator conditions based on anno column name
        observeEvent(input$variable_diag, {
            req(is(get_se(), "NxtSE"))
            req(input$variable_diag != "(none)")
            colData <- colData(get_se())
            req(input$variable_diag %in% colnames(colData))

            if(!is(colData[,input$variable_diag], "factor")) {
                output$warning_diag <- renderText(
                    "Contrast must be performed on discrete categories")
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = "(none)"
                )
            } else {
                output$warning_diag <- renderText("")
                col_levels <- levels(colData[,input$variable_diag])
                updateSelectInput(
                    session = session, inputId = "nom_diag", 
                    choices = c("(none)", col_levels), 
                    selected = "(none)"
                )
                updateSelectInput(
                    session = session, inputId = "denom_diag", 
                    choices = c("(none)", col_levels), 
                    selected = "(none)"
                )
            }
        })

        # Reset to default
        observeEvent(input$clear_diag, {           
            if(is_valid(is(get_se(), "NxtSE"))) {
                colData <- colData(get_se())
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = "(none)"
                )
            } else {
                updateSelectInput(session = session, inputId = "variable_diag",
                    choices = c("(none)"), selected = "(none)")
            }
            
            updateSelectInput(session = session, inputId = "nom_diag", 
                choices = c("(none)"), selected = "(none)")
            updateSelectInput(session = session, inputId = "denom_diag", 
                choices = c("(none)"), selected = "(none)")

            updateSwitchInput(session = session, inputId = "NMD_diag",
                value = FALSE)
        })

        observeEvent(input$clear_selected, {
            req(input$clear_selected)
            settings_Diag$selected <- NULL
        })
    
        get_ggplot <- reactive({
            settings_Diag$ggplot
        })
        spModule <- vis_ggplot_server("scatterSave", get_ggplot, volumes)
    
        return(settings_Diag)        
    })
}

.get_volcano_data_FCunits <- function(res) {
    if("log2FoldChange" %in% colnames(res)) {
        return("log2FoldChange")
    } else if("logFC" %in% colnames(res)) {
        return("logFC")
    } else if("MLE_LFC" %in% colnames(res)) {
        return("MLE_LFC")
    }
}

.get_volcano_data_sigunits <- function(res) {
    if("pvalue" %in% colnames(res)) {
        return(c("pvalue", "padj"))       # DESeq2
    } else if("P.Value" %in% colnames(res)) {
        return(c("P.Value", "adj.P.Val")) # limma or DoubleExpSeq
    } else if("PValue" %in% colnames(res)) {
        return(c("PValue", "FDR"))        # edgeR
    }
}

.get_unified_volcano_data <- function(res) {
    res <- as.data.table(res)
    xunits <- .get_volcano_data_FCunits(res)
    yunits <- .get_volcano_data_sigunits(res)
    df.volc <- data.frame(
        EventName = res$EventName, 
        EventType = res$EventType, 
        NMD_direction = .getNMDcode(res$flags),
        logFC = res[, get(xunits)],
        pvalue = res[, get(yunits[1])],
        FDR = res[, get(yunits[2])]
    )
    return(df.volc)
}

server_vis_volcano <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Volc <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
            randomVar <- runif(1)
            settings_Volc$trigger <- randomVar
        })

        useDE_r <- visFilter_server("volcano", get_de, rows_all, rows_selected)

        observeEvent(rows_selected(), {
            settings_Volc$selected <- rows_selected()
        }, ignoreNULL = FALSE)

        observe({
            req(settings_Volc$trigger)
            req(is(get_se(), "NxtSE"))
            req(get_de())
            
            res <- useDE_r()
            req(res)
            req(nrow(res) > 0)
            req(all(res$EventName %in% rownames(get_se())))

            df.volc <- as.data.frame(res)
            if(input$NMD_volc) {
                df.volc <- df.volc[df.volc$NMD_direction != 0, ]
                df.volc$logFC <- df.volc$logFC * df.volc$NMD_direction
            }
            req(nrow(df.volc) > 0)

            xunits <- .get_volcano_data_FCunits(res)           
            selected <- settings_Volc$selected
            if(is_valid(selected)) {
                df.volc$selected <- 
                    (df.volc$EventName %in% get_de()$EventName[selected])
            } else {
                df.volc$selected <- FALSE
            }

            settings_Volc$plot_ini <- TRUE

            if(input$adjP_volc) {
                p <- ggplot(df.volc, aes(
                        x = get("logFC"), y = -log10(get("FDR")),
                        key = get("EventName"), text = get("EventName"), 
                        colour = get("selected")))           
            } else {
                p <- ggplot(df.volc, aes(
                        x = get("logFC"), y = -log10(get("pvalue")),
                        key = get("EventName"), text = get("EventName"), 
                        colour = get("selected")))               
            }

            p <- p + geom_point() + 
                scale_color_manual(values = c("black", "red"))

            if(input$facet_volc) {
                p <- p + facet_wrap(vars(get("EventType")))
            }
            if(xunits %in% c("log2FoldChange")) {
                formatted_units <- "Log2 Fold Change"
            } else if(xunits %in% c("logFC")) {
                formatted_units <- "Log Fold Change"
            } else {
                formatted_units <- "MLE Log2 Fold Change"
            }
            if(input$NMD_volc) {
                p <- p + labs(x = paste(formatted_units, "NMD substrate"))
            } else {
                p <- p + labs(x = formatted_units)            
            }
            if(input$adjP_volc) {
                p <- p + labs(y = "Adjusted P Value (-log10)")
            } else {
                p <- p + labs(y = "Nominal P Value (-log10)")            
            }
            
            p <- p + labs(color = "Selected")
            settings_Volc$ggplot <- p

            withProgress(message = 'Rendering plot...', value = 0, {
                py <- ggplotly(
                    p, tooltip = "text",
                    source = "plotly_volcano",
                    type = "scatter_gl"
                ) %>% toWebGL() %>% layout(dragmode = "select")

                # Add hoveron entry
                py$x$data <- lapply(py$x$data, function(x) {
                    x$hoveron <- NULL
                    x
                })
                settings_Volc$final_plot <- py

                if(packageVersion("plotly") >= "4.9.0") {
                    plotly::event_register(
                        settings_Volc$final_plot, "plotly_click")
                    plotly::event_register(
                        settings_Volc$final_plot, "plotly_selected")
                }
            })
            
        })

        output$plot_volc <- renderPlotly({
            validate(need(is(get_se(), "NxtSE"), "Load Experiment first"))
            validate(need(get_de(), "Perform DE Analysis first"))

            res <- useDE_r()
            validate(need(nrow(res) > 0, 
                "No events found. Consider relaxing some filters"
            ))
            validate(need(all(res$EventName %in% rownames(get_se())),
                "Filtered NxtSE does not match DE results"
            ))

            df.volc <- as.data.frame(res)
            if(input$NMD_volc) df.volc <- df.volc[df.volc$NMD_direction != 0, ]
            validate(need(nrow(df.volc) > 0, 
                "No events found. Consider relaxing some filters"
            ))
            
            req(settings_Volc$final_plot)
            settings_Volc$final_plot
        })

        # Reactive click
        settings_Volc$plotly_click <- reactive({
            plot_exist <- settings_Volc$plot_ini
            if(plot_exist) 
                event_data("plotly_click", source = "plotly_volcano")
        })
        observeEvent(settings_Volc$plotly_click(), {
            req(settings_Volc$plotly_click())
            click <- settings_Volc$plotly_click()
            click.id <- which(get_de()$EventName == click$key)
            req(click.id)

            selected <- settings_Volc$selected

            if(click.id %in% selected && input$reverse_select) {
                selected <- selected[-which(selected == click.id)]
            } else if(!input$reverse_select) {
                selected <- c(selected, click.id)
            }
            settings_Volc$selected <- selected
        })

        # Reactive brush
        settings_Volc$plotly_brush <- reactive({
            plot_exist <- settings_Volc$plot_ini
            if(plot_exist)
                event_data("plotly_selected", source = "plotly_volcano")
        })
        observeEvent(settings_Volc$plotly_brush(), {
            req(settings_Volc$plotly_brush())
            brush <- settings_Volc$plotly_brush()
            brush.id <- which(get_de()$EventName %in% brush$key)
            req(brush.id)

            selected <- settings_Volc$selected
            if(!input$reverse_select) {
                selected <- union(selected, brush.id)           
            } else {
                selected <- setdiff(selected, brush.id)
            }
            settings_Volc$selected <- selected
        })

        observeEvent(input$clear_volc, {
            updateSwitchInput(session = session, 
                inputId = "facet_volc", value = FALSE)
            updateSwitchInput(session = session, 
                inputId = "adjP_volc", value = TRUE)
            updateSwitchInput(session = session, 
                inputId = "NMD_volc", value = FALSE)
        })

        observeEvent(input$clear_selected, {
            req(input$clear_selected)
            settings_Volc$selected <- NULL
        })

        get_ggplot <- reactive({
            settings_Volc$ggplot
        })
        spModule <- vis_ggplot_server("volcanoSave", get_ggplot, volumes)

        return(settings_Volc)
    })
}

server_vis_heatmap <- function(
        id, refresh_tab, volumes, get_se, get_de, get_go, 
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Heat <- setreactive_Diag()

        useDE_r <- visFilter_server("heatmap", get_de, rows_all, rows_selected)

        # Tab refresh
        # - update conditions
        # - update GO terms
        observeEvent(refresh_tab(), {
            req(refresh_tab())
            randomVar <- runif(1)
            settings_Heat$trigger <- randomVar
            # print(randomVar)
        })
        
        observeEvent(settings_Heat$trigger, {
            req(settings_Heat$trigger)
            req(is(get_se(), "NxtSE"))
            colData <- colData(get_se())
            if(
                    is_valid(input$anno_col_heat) && 
                    all(input$anno_col_heat %in% colnames(colData))
            ) {
                selected <- isolate(input$anno_col_heat)
                updateSelectInput(session = session, 
                    inputId = "anno_col_heat", 
                    choices = colnames(colData(get_se())), 
                    selected = selected)
                if(
                    is_valid(input$anno_col_heat_sort) && 
                    input$anno_col_heat_sort %in% colnames(colData)
                ) {
                    selected2 <- isolate(input$anno_col_heat_sort)
                    updateSelectInput(session = session, 
                        inputId = "anno_col_heat_sort", 
                        choices = c("(none)", selected), 
                        selected = selected2)
                } else {
                    updateSelectInput(session = session, 
                        inputId = "anno_col_heat_sort", 
                        choices = c("(none)", selected), 
                        selected = "(none)")
                }
            } else {
                updateSelectInput(session = session, 
                    inputId = "anno_col_heat", 
                    choices = colnames(colData(get_se())), 
                    selected = NULL)
                updateSelectInput(session = session, 
                    inputId = "anno_col_heat_sort", 
                    choices = c("(none)"), 
                    selected = "(none)")
            }

            # Update annotation column names in selection
            req(get_go())
            goTerms <- as.character(get_go()$Term)
            if(length(goTerms) > 50) goTerms <- goTerms[seq_len(50)]
            if(
                    is_valid(input$GO_heat) && 
                    input$GO_heat %in% goTerms
            ) {
                selectedOption <- isolate(input$GO_heat)
                updateSelectInput(
                    session = session, inputId = "GO_heat", 
                    choices = c("(none)", goTerms), 
                    selected = selectedOption
                )
            } else {
                updateSelectInput(
                    session = session, inputId = "GO_heat", 
                    choices = c("(none)", goTerms)
                )
            }
        })
        
        # Enable sample annotations
        observeEvent(input$anno_col_heat, {
            req(is(get_se(), "NxtSE"))
            colData <- colData(get_se())
            if(
                    is_valid(input$anno_col_heat) && 
                    all(input$anno_col_heat %in% colnames(colData))
            ) {
                selected <- isolate(input$anno_col_heat)
                if(
                    is_valid(input$anno_col_heat_sort) && 
                    input$anno_col_heat_sort %in% colnames(colData)
                ) {
                    selected2 <- isolate(input$anno_col_heat_sort)
                    updateSelectInput(session = session, 
                        inputId = "anno_col_heat_sort", 
                        choices = c("(none)", selected), 
                        selected = selected2)
                } else {
                    updateSelectInput(session = session, 
                        inputId = "anno_col_heat_sort", 
                        choices = c("(none)", selected), 
                        selected = "(none)")
                }
            } else {
                updateSelectInput(session = session, 
                    inputId = "anno_col_heat_sort", 
                    choices = c("(none)"), 
                    selected = "(none)")
            }
        })

        observe({
            res <- get_de()
            if(!is_valid(res)) {
                settings_Heat$eventsGO <- NULL
            } else if(is_valid(input$GO_heat)) {
                goInfo <- get_go()
                go_id <- goInfo$go_id[match(input$GO_heat, goInfo$Term)]
                settings_Heat$eventsGO <- subset_EventNames_by_GO(
                    res$EventName, go_id, get_se())
            } else {
                settings_Heat$eventsGO <- res$EventName
            }
        })
        
        observe({
            req(settings_Heat$trigger)
            req(is(get_se(), "NxtSE"))
            req(useDE_r())

            res <- useDE_r()
            req(nrow(res) > 0)
            res <- res[, c("EventName", "EventType"), with = FALSE]
            res <- res[get("EventName") %in% settings_Heat$eventsGO]
            req(nrow(res) > 0)            
            req(all(res$EventName %in% rownames(get_se())))
            
            colData <- as.data.frame(colData(get_se()))
            
            rowLim <- input$slider_num_events_heat
            if(nrow(res) > rowLim) {
                for(rNum in seq(rowLim, nrow(res), by = rowLim)) {
                    resSubset <- res[seq_len(rNum)]
                    if(input$mode_heat == "PSI") {
                        mat <- makeMatrix(get_se(), resSubset$EventName,
                            rownames(colData), "PSI")
                    } else if(input$mode_heat == "Logit") {
                        mat <- makeMatrix(get_se(), resSubset$EventName,
                            rownames(colData), "logit")
                    } else {
                        mat <- makeMatrix(get_se(), resSubset$EventName,
                            rownames(colData), "Z-score")
                    }
                    if(nrow(mat) >= rowLim) {
                        mat <- mat[seq_len(rowLim),]
                        break
                    }
                }            
            } else {
                if(input$mode_heat == "PSI") {
                    mat <- makeMatrix(get_se(), res$EventName,
                        rownames(colData), "PSI")
                } else if(input$mode_heat == "Logit") {
                    mat <- makeMatrix(get_se(), res$EventName,
                        rownames(colData), "logit")
                } else {
                    mat <- makeMatrix(get_se(), res$EventName,
                        rownames(colData), "Z-score")
                }
            }
            settings_Heat$mat <- mat
        })

        observe({
            req(get_se())
            req(settings_Heat$mat)
            req(nrow(settings_Heat$mat) > 0)

            mat <- settings_Heat$mat
            
            colData <- as.data.frame(colData(get_se()))

            colors.df <- RColorBrewer::brewer.pal.info
            color.index <- which(rownames(colors.df) == input$color_heat)
            color <- grDevices::colorRampPalette(
                rev(RColorBrewer::brewer.pal(
                    n = colors.df$maxcolors[color.index],
                    name = rownames(colors.df)[color.index])
                )
            )
            color_vec <- color(100)
            # Hopefully the fixed filtering in limma pipeline will also fix the 
            #   NA issues here:
            na.exclude <- (rowSums(!is.na(mat)) == 0)
            if(any(na.exclude == TRUE)) {
                output$warning_heat <- renderText({
                    "Some events have been excluded due to NA values"
                    # paste(rownames(mat)[which(na.exclude)])
                })
                mat <- mat[-which(na.exclude),]
            }

            if(
                    is_valid(input$anno_col_heat) && 
                    all(input$anno_col_heat %in% colnames(colData))
            ) {
                if(
                    is_valid(input$anno_col_heat_sort) &&
                    input$anno_col_heat_sort %in% colnames(colData)
                ) {
                    new_order <- order(
                        colData[, input$anno_col_heat_sort],
                        decreasing = input$anno_col_heat_sort_order
                    )
                    mat <- mat[, new_order]
                    colData_sorted <- colData[new_order, ]
                    
                    settings_Heat$ggplot <- pheatmap(
                        mat, color = color_vec, 
                        annotation_col = colData_sorted[, 
                            input$anno_col_heat, drop=FALSE],
                        cluster_cols = FALSE
                    )
                    # settings_Heat$ggplot <- heatmaply::ggheatmap(
                        # mat, color = color, 
                        # col_side_colors = colData_sorted[, 
                            # input$anno_col_heat, drop=FALSE],
                        # dendrogram = "row"
                    # )
                    settings_Heat$final_plot <- heatmaply::heatmaply(
                        mat, color = color, 
                        col_side_colors = colData_sorted[, 
                            input$anno_col_heat, drop=FALSE],
                        dendrogram = "row"
                    )
                } else {
                    settings_Heat$ggplot <- pheatmap(
                        mat, color = color_vec, 
                        annotation_col = colData[, 
                            input$anno_col_heat, drop=FALSE]
                    )
                    # settings_Heat$ggplot <- heatmaply::ggheatmap(
                        # mat, color = color, 
                        # col_side_colors = colData[, 
                            # input$anno_col_heat, drop=FALSE]
                    # )
                    settings_Heat$final_plot <- heatmaply::heatmaply(
                        mat, color = color, 
                        col_side_colors = colData[, 
                            input$anno_col_heat, drop=FALSE]
                    )
                }

            } else {
                settings_Heat$ggplot <- pheatmap(
                    mat, color = color_vec
                )
                # settings_Heat$ggplot <- heatmaply::ggheatmap(
                    # mat, color = color)
                settings_Heat$final_plot <- heatmaply::heatmaply(
                    mat, color = color)
            }      
            
        })

        output$plot_heat <- renderPlotly({
            validate(need(is(get_se(), "NxtSE"), "Load Experiment first"))
            validate(need(get_de(), "Perform DE analysis first"))

            res <- useDE_r()
            validate(need(nrow(res) > 0,
                "No events to display. Suggest relaxing some filters"
            ))
            res <- res[, c("EventName", "EventType"), with = FALSE]
            res <- res[get("EventName") %in% settings_Heat$eventsGO]
            validate(need(nrow(res) > 0,
                "No events to display after gene ontology filtering"
            ))

            validate(need(all(res$EventName %in% rownames(get_se())),
                "Filtered NxtSE does not match DE results"
            ))

            validate(need(settings_Heat$mat, 
                "No events with finite values to display"
            ))
            validate(need(nrow(settings_Heat$mat) > 0, 
                "No events with finite values to display"
            ))

            req(settings_Heat$final_plot)
            settings_Heat$final_plot
        })

        get_ggplot <- reactive({
            settings_Heat$ggplot
        })
        spModule <- vis_ggplot_server("heatSave", get_ggplot, volumes)
        
    })
}

.getNMDcode <- function(vec) {
    IncIsNMD <- grepl("Inc-NMD", vec)
    ExcIsNMD <- grepl("Exc-NMD", vec)
    return(ifelse(IncIsNMD, 1, ifelse(ExcIsNMD, -1, 0)))
}