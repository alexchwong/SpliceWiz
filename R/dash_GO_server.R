server_GO <- function(
        id, refresh_tab, nxtse_path, get_se, get_de,
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
       
        # Plot GO
        observeEvent(input$plot_GO, {
            
            output$plot_GO <- renderPlotly({
                validate(need(get_se(), "Load Experiment first"))
                validate(need(get_de(), "Load DE Analysis first"))

                ref_path <- file.path(settings_GO$nxtse_path, "Reference")
                ontFile <- file.path(ref_path, "fst/Ontology.fst")
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
                res <- .get_unified_volcano_data(res_ET)
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
                    res <- res[get("LogFC") > 0]
                } else if(input$direction_GO == "Down") {
                    res <- res[get("LogFC") < 0]
                }
                
                validate(need(nrow(res) > 0, "Zero differential events"))
                selectedEvents <- res$EventType
                
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
               
                res <- goASE(
                    selectedEvents, universeEvents,
                    ref_path, ontologyType
                )
                
                print(res)
            })
        })
        
    })
    
}