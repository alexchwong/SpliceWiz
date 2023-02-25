server_GO <- function(
        id, refresh_tab, get_se, get_de, volumes,
        rows_all, rows_selected
) {

    moduleServer(id, function(input, output, session) {
        settings_GO <- setreactive_GO()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
        })

        observe({
            shinyFileSave(input, "GO_export_geneId", 
                roots = volumes(), session = session, filetypes = c("txt"))    
        })
        observeEvent(input$GO_export_geneId, {
            selectedfile <- parseSavePath(volumes(), input$GO_export_geneId)
            req(selectedfile$datapath)
            
            # Save gene id's to file
            req(settings_GO$gene_ids)
            fwrite(list(settings_GO$gene_ids), selectedfile$datapath)
        })

        observe({
            shinyFileSave(input, "GO_export_univId", 
                roots = volumes(), session = session, filetypes = c("txt"))    
        })
        observeEvent(input$GO_export_univId, {
            selectedfile <- parseSavePath(volumes(), input$GO_export_univId)
            req(selectedfile$datapath)
            
            # Save (bkgd) gene id's to file
            req(settings_GO$univ_ids)
            fwrite(list(settings_GO$univ_ids), selectedfile$datapath)
        })

        observeEvent(input$perform_GO, {
            output$warning_GO <- renderText(isolate({
                validate(need(is(get_se(), "NxtSE"), "Load NxtSE first"))
                validate(need(get_de(), "Load DE Analysis first"))

                validate(need(ref(get_se())$ontology,
                    "No gene ontology found for this NxtSE object"))
                
                # Get volcano data
                
                res_bkgd <- as.data.table(.get_unified_volcano_data(get_de()))
                xunits <- .get_volcano_data_FCunits(res_bkgd)
                res_bkgd$All <- seq_len(nrow(res_bkgd)) %in% rows_all()
                res_bkgd$Selected <- seq_len(nrow(res_bkgd)) %in% rows_selected()

                res_all <- res_bkgd[get("All") == TRUE]
               if(is_valid(input$EventType_GO)) {
                    res_all <- res_all[get("EventType") %in% input$EventType_GO]
                } else {
                    res_all <- res_all
                }
                
                validate(need(nrow(res_all) > 0, "Zero differential events"))
                # Filter for Top N events or significant events
                if(input$threshType_GO == "Top events by p-value") {
                    req(input$topN_GO)
                    res <- res_all[seq_len(input$topN_GO)]
                } else if(input$threshType_GO == "Nominal P value") {
                    req(input$pvalT_GO)
                    res <- res_all[get("pvalue") <= input$pvalT_GO]
                } else if(input$threshType_GO == "Adjusted P value") {
                    req(input$pvalT_GO)
                    res <- res_all[get("FDR") <= input$pvalT_GO]                
                } else if(input$threshType_GO == "Highlighted events") {
                    res <- res_all[get("Selected") == TRUE]
                }

                validate(need(nrow(res) > 0, "Zero differential events"))
                if(input$direction_GO == "Up") {
                    res <- res[get(xunits) > 0]
                } else if(input$direction_GO == "Down") {
                    res <- res[get(xunits) < 0]
                }
                
                validate(need(nrow(res) > 0, "Zero differential events"))
                selectedEvents <- res$EventName
                
                # debug
                # print(selectedEvents)
                
                # Get Universe
                universeEvents <- res_all$EventName
                if(input$universe_GO == "All Genes") {
                    # a signal to use all genes instead            
                    universeEvents <- NULL 
                }

                ontologyType <- "BP"
                if(input$category_GO == "Molecular Function") {
                    ontologyType <- "MF"
                } else if(input$category_GO == "Cellular Compartment") {
                    ontologyType <- "CC"
                }
                
                withProgress(message = 'Performing GO analysis...', value = 0, {
                    # Store filtered volcano-data
                    settings_GO$filteredVolc <- res
                    
                    # Get gene_ids for ASEs (for optional save to file)
                    geneIds <- extract_gene_ids_for_GO(
                        selectedEvents,
                        universeEvents,
                        get_se()
                    )
                    settings_GO$gene_ids <- geneIds$genes
                    settings_GO$univ_ids <- geneIds$universe
                    
                    # Generate GO
                    settings_GO$resGO <- .format_GO_result(
                        .ora_internal(
                            ref(get_se())[["ontology"]], 
                            geneIds$genes, geneIds$universe,
                            ontologyType, pAdjustMethod = "BH"
                        )
                    )
                })
                
                p <- .generate_ggplot_GO(settings_GO$resGO)
                settings_GO$final_plot <- ggplotly(
                    p, tooltip = "text"
                )
                
                "GO analysis complete"
            }))
        })

        output$plot_GO <- renderPlotly({
            validate(need(nrow(settings_GO$resGO) > 0, "Zero results"))
            print(settings_GO$final_plot)
        })
        
        return(settings_GO)
    })
    
}

