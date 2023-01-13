# Order of DE filtration
# get_de()
# - get_de()[rows_all]
# - get_de()[rows_all][seq_len(Top N rows)]

# New filtration approach
# get_de()
# - get_de()[rows_all]
#   - {option to filter by padj, pvalue, or top n rows}

server_vis_diag <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Diag <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
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

        # Reactive to generate filtered DE object
        observe({
            req(get_de())
            tmpres <- as.data.table(
                .get_unified_volcano_data(get_de()[rows_all(),]))
            if(input$filterType_diag == "Adjusted P value") {
                settings_Diag$useDE <- tmpres[get("FDR") <= input$pvalT_diag]
            } else if(input$filterType_diag == "Nominal P value") {
                settings_Diag$useDE <- tmpres[get("pvalue") <= input$pvalT_diag]
            } else if(input$filterType_diag == "Top N results") {
                if(input$topN_diag < nrow(settings_Diag$useDE)) {
                    settings_Diag$useDE <- tmpres[seq_len(input$topN_diag)]
                } else {
                    settings_Diag$useDE <- tmpres
                }
            }
        })

        # Update local rows_selected with that of global
        observeEvent(rows_selected(), {
            settings_Diag$selected <- rows_selected()
        }, ignoreNULL = FALSE)
    
        output$plot_diag <- renderPlotly({
            validate(need(is(get_se(), "NxtSE"), "Load Experiment first"))
            validate(need(settings_Diag$useDE, "Perform DE Analysis first"))
            validate(need(input$variable_diag, 
                "Select conditions and contrasts"))
            validate(need(input$nom_diag, 
                "Select conditions and contrasts"))
            validate(need(input$denom_diag, 
                "Select conditions and contrasts"))
            validate(need(input$variable_diag != "(none)", 
                "Select conditions and contrasts"))
            validate(need(input$nom_diag != "(none)", 
                "Select conditions and contrasts"))
            validate(need(input$denom_diag != "(none)", 
                "Select conditions and contrasts"))

            # Filter DE by EventType; fetch diag object
            res <- settings_Diag$useDE
            if(is_valid(input$EventType_diag)) {
                res <- res[get("EventType") %in% input$EventType_diag]
            }
            withProgress(message = 'Calculating mean PSIs...', value = 0, {
                df.diag <- makeMeanPSI(
                    get_se(), res$EventName, input$variable_diag, 
                    list(input$nom_diag, input$denom_diag)
                )
                colnames(df.diag)[seq(2,3)] <- c("nom", "denom")
            })
            
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
            
            # Generate ggplot object
            settings_Diag$plot_ini <- TRUE
            if(input$NMD_diag == TRUE) {
                df.diag             <- df.diag[df.diag$NMD_direction != 0, ]
                df.diag$nom_NMD     <- ifelse(df.diag$NMD_direction == 1, 
                                        df.diag$nom, df.diag$denom)
                df.diag$denom_NMD   <- ifelse(df.diag$NMD_direction == -1, 
                                        df.diag$nom, df.diag$denom)
                                        
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
            
            # Record ggplot / plotly objects into settings_Diag
            settings_Diag$ggplot <- p
            settings_Diag$final_plot <- ggplotly(
                p, tooltip = "text",
                source = "plotly_diagonal"
            ) %>% layout(
                dragmode = "lasso",
                yaxis = list(scaleanchor="x", scaleratio=1)
            )
            if(packageVersion("plotly") >= "4.9.0") {
                plotly::event_register(
                    settings_Diag$final_plot, "plotly_click")
                plotly::event_register(
                    settings_Diag$final_plot, "plotly_selected")
            }
            
            withProgress(message = 'Rendering plot...', value = 0, {
                print(settings_Diag$final_plot)
            })
        })

        # Output ggplot to RStudio plot window
        observeEvent(input$output_plot_diag, {
            req(settings_Diag$ggplot)
            print(isolate(settings_Diag$ggplot))
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

            if(click.id %in% selected) {
                selected <- selected[-which(selected == click.id)]
            } else {
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
            selected <- unique(c(selected, brush.id))
            settings_Diag$selected <- selected
        })
    
        # Update nominator / denominator conditions based on anno column name
        observeEvent(input$variable_diag, {
            req(get_se())
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
            updateSelectInput(session = session, 
                "EventType_diag", selected = NULL)
            updateSelectInput(session = session, 
                "filterType_diag", selected = "Adjusted P value")
            shinyWidgets::updateSliderTextInput(session = session, 
                "topN_diag", selected = 500)
            shinyWidgets::updateSliderTextInput(session = session, 
                "pvalT_diag", selected = 0.05)
            
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
        })
    
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
        })

        # Reactive to generate filtered DE object
        observe({
            req(get_de())
            tmpres <- as.data.table(
                .get_unified_volcano_data(get_de()[rows_all(),]))
            if(input$filterType_volc == "Adjusted P value") {
                settings_Volc$useDE <- tmpres[get("FDR") <= input$pvalT_volc]
            } else if(input$filterType_volc == "Nominal P value") {
                settings_Volc$useDE <- tmpres[get("pvalue") <= input$pvalT_volc]
            } else if(input$filterType_volc == "Top N results") {
                if(input$topN_volc < nrow(settings_Volc$useDE)) {
                    settings_Volc$useDE <- tmpres[seq_len(input$topN_volc)]
                } else {
                    settings_Volc$useDE <- tmpres
                }
            }
        })

        observeEvent(rows_selected(), {
            settings_Volc$selected <- rows_selected()
        }, ignoreNULL = FALSE)

        output$plot_volc <- renderPlotly({
            validate(need(is(get_se(), "NxtSE"), "Load Experiment first"))
            validate(need(settings_Volc$useDE, "Perform DE Analysis first"))

            selected <- settings_Volc$selected

            res <- settings_Volc$useDE
            if(is_valid(input$EventType_volc)) {
                res <- res[get("EventType") %in% input$EventType_volc]
            }

            xunits <- .get_volcano_data_FCunits(get_de())
            df.volc <- as.data.frame(res)
            
            if(is_valid(selected)) {
                df.volc$selected <- 
                    (df.volc$EventName %in% get_de()$EventName[selected])
            } else {
                df.volc$selected <- FALSE
            }
            if(input$NMD_volc) {
                df.volc <- df.volc[df.volc$NMD_direction != 0, ]
                df.volc$logFC <- df.volc$logFC * df.volc$NMD_direction
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
            settings_Volc$final_plot <- ggplotly(
                p, tooltip = "text",
                source = "plotly_volcano"
            ) %>% layout(dragmode = "select")
            if(packageVersion("plotly") >= "4.9.0") {
                plotly::event_register(
                    settings_Volc$final_plot, "plotly_click")
                plotly::event_register(
                    settings_Volc$final_plot, "plotly_selected")
            }
            withProgress(message = 'Rendering plot...', value = 0, {
                print(settings_Volc$final_plot)
            })
        })

        observeEvent(input$output_plot_volc, {
            req(settings_Volc$ggplot)
            print(isolate(settings_Volc$ggplot))
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

            if(click.id %in% selected) {
                selected <- selected[-which(selected == click.id)]
            } else {
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
            selected <- unique(c(selected, brush.id))
            settings_Volc$selected <- selected
        })

        observeEvent(input$clear_volc, {
            updateSelectInput(session = session, "EventType_volc", 
                selected = NULL)
            shinyWidgets::updateSliderTextInput(session = session, 
                "number_events_volc", selected = 1000)
        })
        
        return(settings_Volc)
    })
}

server_vis_heatmap <- function(
        id, refresh_tab, volumes, get_se, get_de, get_go, nxtse_path,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Heat <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
            req(get_se())
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
            goTerms <- get_go()$Term
            if(
                    is_valid(input$GO_heat) && 
                    input$GO_heat %in% goTerms
            ) {
                selectedOption <- isolate(input$GO_heat)
                updateSelectInput(
                    session = session, inputId = "GO_heat", 
                    choices = goTerms, 
                    selected = selectedOption
                )
            } else {
                updateSelectInput(
                    session = session, inputId = "GO_heat", 
                    choices = goTerms
                )
            }
        })
        observeEvent(input$anno_col_heat, {
            req(get_se())
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
        
        # Reactive to generate filtered DE object
        observe({
            req(get_de())
            tmpres <- as.data.table(
                .get_unified_volcano_data(get_de()[rows_all(),]))
            if(input$filterType_heat == "Adjusted P value") {
                tmpres2 <- tmpres[get("FDR") <= input$pvalT_heat]
            } else if(input$filterType_heat == "Nominal P value") {
                tmpres2 <- tmpres[get("pvalue") <= input$pvalT_heat]
            } else if(input$filterType_heat == "Top N results") {
                if(input$topN_heat < nrow(settings_Heat$useDE)) {
                    tmpres2 <- tmpres[seq_len(input$topN_heat)]
                } else {
                    tmpres2 <- tmpres
                }
            }
            settings_Heat$useDE <- tmpres2[, c("EventType"), with = FALSE]
        })

        output$plot_heat <- renderPlotly({
            
            validate(need(is(get_se(), "NxtSE"), "Load Experiment first"))
            validate(need(settings_Heat$useDE, "Perform DE Analysis first"))

            res <- settings_Heat$useDE
            
            # Filter by highlighted events or GO category
            if(input$secondFilter_heat == "Highlighted (selected) events") {
                selected <- rows_selected()
                res <- res[get("EventName") %in% get_de()$EventName[selected]]
            } else if(input$secondFilter_heat == 
                "Top Gene Ontology Categories")
            {
                # filter by selected GO category
                validate(need(get_go(), 
                    "Run Gene Ontology analysis first"))
                goInfo <- get_go()[get("Term") == input$GO_heat]
                go_id <- goInfo$go_id
                events <- .subset_EventNames_by_GO(res$EventName, go_id,
                    nxtse_path())
                res <- res[get("EventName") %in% events]
            } else {
                # do nothing
            }
            
            if(nrow(res) > input$slider_num_events_heat) {
                res <- res[seq_len(input$slider_num_events_heat)]
            }
            
            validate(need(nrow(res) > 0, "No events to plot"))

            colData <- as.data.frame(colData(get_se()))

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

            validate(need(nrow(mat) > 0 & ncol(mat) > 0, 
                "No Events with sufficient finite PSI values to draw heatmap"))

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
                settings_Heat$final_plot <- heatmaply::heatmaply(
                    mat, color = color)
            }      
            print(settings_Heat$final_plot)
        })

    })
}

.getNMDcode <- function(vec) {
    IncIsNMD <- grepl("Inc-NMD", vec)
    ExcIsNMD <- grepl("Exc-NMD", vec)
    return(ifelse(IncIsNMD, 1, ifelse(ExcIsNMD, -1, 0)))
}