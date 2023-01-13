server_GO <- function(
        id, refresh_tab, nxtse_path, get_de, volumes,
        rows_all, rows_selected
) {

    moduleServer(id, function(input, output, session) {
        settings_GO <- setreactive_GO()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
        })

        observeEvent(nxtse_path(), {
            req(nxtse_path())
            settings_GO$nxtse_path <- nxtse_path()
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
                validate(need(settings_GO$nxtse_path, "Load NxtSE first"))
                validate(need(get_de(), "Load DE Analysis first"))

                validate(need(dir.exists(settings_GO$nxtse_path),
                    "NxtSE directory does not exist"
                ))
                
                ontFile <- file.path(reference_path, "fst/Ontology.fst")
                validate(need(file.exists(ontFile),
                    paste(ontFile, "GO reference does not exist")
                ))
                # Get volcano data
                res_all <- as.data.table(get_de()[rows_all(),])
                if(is_valid(input$EventType_GO)) {
                    res_ET <- res_all[get("EventType") %in% input$EventType_GO]
                } else {
                    res_ET <- res_all
                }
                res <- as.data.table(.get_unified_volcano_data(res_ET))
                xunits <- .get_volcano_data_FCunits(res)

                validate(need(nrow(res) > 0, "Zero differential events"))
                # Filter for Top N events or significant events
                if(input$threshType_GO == "Top N results") {
                    req(input$topN_GO)
                    res <- res[seq_len(input$topN_GO)]
                } else if(input$threshType_GO == "Nominal P value") {
                    req(input$pvalT_GO)
                    res <- res[get("pvalue") <= input$pvalT_GO]
                } else {
                    req(input$pvalT_GO)
                    res <- res[get("FDR") <= input$pvalT_GO]                
                }

                validate(need(nrow(res) > 0, "Zero differential events"))
                if(input$direction_GO == "Up") {
                    res <- res[get(xunits) > 0]
                } else if(input$direction_GO == "Down") {
                    res <- res[get(xunits) < 0]
                }
                
                validate(need(nrow(res) > 0, "Zero differential events"))
                selectedEvents <- res$EventName
                
                # Get Universe
                universeEvents <- res_all$EventName
                if(input$universe_GO == "Selected ASE Modality") {
                    universeEvents <- res_ET$EventName
                } else if(input$universe_GO == "All Genes") {
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
                    geneIds <- .extract_gene_ids_for_GO(
                        selectedEvents,
                        universeEvents,
                        settings_GO$nxtse_path
                    )
                    settings_GO$gene_ids <- geneIds$genes
                    settings_GO$univ_ids <- geneIds$universe
                    
                    # Generate GO
                    settings_GO$resGO <- .format_GO_result(
                        .ora_internal(
                            settings_GO$nxtse_path, 
                            geneIds$genes, geneIds$universe,
                            ontologyType, pAdjustMethod = "BH"
                        )
                    )
                })
                
                settings_GO$final_plot <- .generate_plotly_GO(settings_GO$resGO)
                
                "GO analysis complete"
            }))
        })

        output$plot_GO <- renderPlotly({
            print(settings_GO$final_plot)
        })
        
        return(settings_GO)
    })
    
}

.generate_plotly_GO_labels <- function(axis_term) {
    if(axis_term == "log10FDR") {
        return("Enrichment FDR (-log10)")
    } else if(axis_term == "nGenes") {
        return("Number of Enriched Genes")
    } else {
        return("Fold Enrichment")
    }
}

.format_GO_result <- function(res) {
    res$go_term <- substr(res$go_term, 1, trim_go_term + 1)
    res[nchar(get("go_term")) > 50, c("go_term") :=
        paste0(substr(get("go_term"), 1, trim_go_term-3), "...")]
    res$Term <- paste(res$go_term, res$go_id, sep = "~")
    res$Term <- factor(res$Term, res$Term, ordered = TRUE)
    
    res$FDR = res$padj
    res$log10FDR = -log10(res$FDR)
    res$nGenes = res$overlap

    return(res)
}

.generate_plotly_GO <- function(
    res,
    plot_x = c("log10FDR", "foldEnrichment", "nGenes"),
    plot_size = c("nGenes", "foldEnrichment", "log10FDR"),
    plot_color = c("foldEnrichment", "nGenes", "log10FDR"),
    filter_n_terms = 20,
    filter_padj = 0.05,
    filter_pvalue = 0.05,
    trim_go_term = 50
) {
    plot_x <- match.arg(plot_x)
    plot_size <- match.arg(plot_size)
    plot_color <- match.arg(plot_color)

    res_use <- res[seq_len(filter_n_terms)]
    res_use <- res_use[get("pval") <= filter_pvalue]    
    res_use <- res_use[get("padj") <= filter_padj]
        
    p <- ggplot(res_use, aes(text = get("Term"))) + 
        geom_segment(data = res_use, mapping = aes(
            x = 0, xend = get(plot_x), 
            y = get("Term"), yend = get("Term"), 
            color = get(plot_color)
        )) +
        geom_point(data = res_use, mapping = aes(
            x = get(plot_x), y = get("Term"), 
            size = get(plot_size), color = get(plot_color)
        )) +
        scale_colour_gradient(low = "blue", high = "red") +
        scale_y_discrete(limits=rev) +
        labs(
            x = .generate_plotly_GO_labels(plot_x),
            y = "Gene Ontology Term",
            color = .generate_plotly_GO_labels(plot_color), 
            size = .generate_plotly_GO_labels(plot_size)
        )
    
    return(ggplotly(
        p, tooltip = "text"
    ))
}