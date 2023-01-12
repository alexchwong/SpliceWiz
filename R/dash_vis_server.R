server_vis_diag <- function(
        id, refresh_tab, volumes, get_se, get_de,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {
        settings_Diag <- setreactive_Diag()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
            output$warning_diag <- renderText({
                validate(need(get_se(), 
                "Please load Differential Expression via 'Analysis' tab"))
                
                "Differential Expression Loaded"
            })
            req(get_se())
            colData <- colData(get_se())
            if(
                    is_valid(input$variable_diag) && 
                    input$variable_diag %in% colnames(colData)
            ) {
                selected <- isolate(input$variable_diag)
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = selected
                )
            } else {
                updateSelectInput(
                    session = session, inputId = "variable_diag", 
                    choices = c("(none)", colnames(colData)), 
                    selected = "(none)"
                )
            }

        })
        observeEvent(rows_selected(), {
            settings_Diag$selected <- rows_selected()
        }, ignoreNULL = FALSE)
    
        output$plot_diag <- renderPlotly({
            # settings_Diag$plot_ini = FALSE
            validate(need(get_se(), "Load Experiment first"))
            validate(need(get_de(), "Load DE Analysis first"))
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

            selected <- settings_Diag$selected      

            num_events <- input$number_events_diag
            res <- as.data.table(get_de()[rows_all(),])
            if(is_valid(input$EventType_diag)) {
                res <- res[get("EventType") %in% input$EventType_diag]
            }
            if(num_events < nrow(res)) {
                res <- res[seq_len(num_events)]
            }
            df.diag <- makeMeanPSI(
                get_se(), res$EventName, input$variable_diag, 
                list(input$nom_diag, input$denom_diag)
            )
            colnames(df.diag)[seq(2,3)] <- c("nom", "denom")
            if(is_valid(settings_Diag$selected)) {
                df.diag$selected <- 
                    (df.diag$EventName %in% get_de()$EventName[selected])
            } else {
                df.diag$selected <- FALSE
            }
            df.diag$NMD_direction <- .getNMDcode(get_de()$flags[
                match(df.diag$EventName, get_de()$EventName)])
            
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
            p <- p + labs(color = "Selected")
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
            print(settings_Diag$final_plot)
        })

        observeEvent(input$output_plot_diag, {
            req(settings_Diag$ggplot)
            print(settings_Diag$ggplot)
        })
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
            # DT::dataTableProxy("DT_DE") %>% DT::selectRows(selected)
        })

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
            # DT::dataTableProxy("DT_DE") %>% DT::selectRows(selected)
        })
    
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

        observeEvent(input$clear_diag, {
            updateSelectInput(session = session, 
                "EventType_diag", selected = NULL)
            shinyWidgets::updateSliderTextInput(session = session, 
                "number_events_diag", selected = 1000)
            
            if(is_valid(get_se())) {
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
        observeEvent(rows_selected(), {
            settings_Volc$selected <- rows_selected()
        }, ignoreNULL = FALSE)

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


        output$plot_volc <- renderPlotly({
            validate(need(get_se(), "Load Experiment first"))
            validate(need(get_de(), "Load DE Analysis first"))

            selected <- settings_Volc$selected

            num_events <- input$number_events_volc
            res <- as.data.table(get_de()[rows_all(),])
            if(is_valid(input$EventType_volc)) {
                res <- res[get("EventType") %in% input$EventType_volc]
            }
            if(num_events < nrow(res)) {
                res <- res[seq_len(num_events)]
            }

            df.volc <- .get_unified_volcano_data(res)
            xunits <- .get_volcano_data_FCunits(res)

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
            print(settings_Volc$final_plot)
        })

        observeEvent(input$output_plot_volc, {
            req(settings_Volc$ggplot)
            print(settings_Volc$ggplot)
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
        id, refresh_tab, volumes, get_se, get_de,
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
        
        output$plot_heat <- renderPlotly({
            
            validate(need(get_se(), "Load Experiment first"))
            validate(need(get_de(), "Load DE Analysis first"))

            if(input$select_events_heat == "Selected") {
                selected <- rows_selected()
            } else if(input$select_events_heat == "Filtered") {
                selected <- rows_all()
                if(length(selected) > input$slider_num_events_heat) {
                    selected <- selected[seq_len(input$slider_num_events_heat)]
                }
            } else {
                selected <- seq_len(min(input$slider_num_events_heat, 
                    nrow(get_de())))
            }

            validate(need(length(selected) > 0, "Select some Events first"))

            colData <- as.data.frame(colData(get_se()))

            if(input$mode_heat == "PSI") {
                mat <- makeMatrix(get_se(), get_de()$EventName[selected],
                rownames(colData), "PSI")
            } else if(input$mode_heat == "Logit") {
                mat <- makeMatrix(get_se(), get_de()$EventName[selected],
                rownames(colData), "logit")
            } else {
                mat <- makeMatrix(get_se(), get_de()$EventName[selected],
                rownames(colData), "Z-score")
            }

            validate(need(nrow(mat) > 0 & ncol(mat) > 0, 
                "No data after filtering results"))

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