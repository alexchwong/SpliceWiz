server_GO <- function(
        id, refresh_tab, get_se, get_de, volumes,
        rows_all, rows_selected
) {

    moduleServer(id, function(input, output, session) {
        settings_GO <- setreactive_GO()

        observeEvent(refresh_tab(), {
            req(refresh_tab())
            settings_GO$trigger <- runif(1)
        })
        
        useDE_r <- visFilter_server("GOfilters", get_de, rows_all, rows_selected)

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

       
        # observe to calculate filtered volcano results
        observe({
            if(!is(get_se(), "NxtSE")) {
                settings_GO$errorMsg <- 
                    "Load NxtSE object first"
            }
            req(is(get_se(), "NxtSE"))
            if(!is_valid(get_de())) {
                settings_GO$errorMsg <- 
                    "Perform differential analysis first"
            }
            req(get_de())
            if(!("ontology" %in% names(ref(get_se())))) {
                settings_GO$errorMsg <- 
                    "No gene ontology found for this NxtSE object"
            }
            req("ontology" %in% names(ref(get_se())))
            req(get_de())
            res_all <- as.data.table(.get_unified_volcano_data(get_de()))
            res <- useDE_r()
            if(!is_valid(res) || nrow(res) <= 0) {
                settings_GO$errorMsg <- "Zero differential events"
            }
            req(res)
            req(nrow(res) > 0)
            xunits <- .get_volcano_data_FCunits(res)
            if(input$direction_GO == "Up") {
                res <- res[get(xunits) > 0]
            } else if(input$direction_GO == "Down") {
                res <- res[get(xunits) < 0]
            }
            if(nrow(res) <= 0) {
                settings_GO$errorMsg <- "Zero differential events"
            }
            
            if(
                !is_valid(res) ||
                nrow(res) == 0           
            ) {
                settings_GO$gene_ids <- NULL
                settings_GO$univ_ids <- NULL
                settings_GO$resGO <- NULL
                settings_GO$ggplot <- NULL
                settings_GO$final_plot <- NULL
            }
            req(nrow(res) > 0)

            if(!all(res$EventName %in% rownames(get_se()))) {
                settings_GO$errorMsg <- 
                    "Filtered NxtSE does not match DE results"
            }
            req(all(res$EventName %in% rownames(get_se())))
            
            settings_GO$filteredVolc <- res

            selectedEvents <- res$EventName

            # Get Universe
            universeEvents <- res_all$EventName
            if(input$universe_GO == "All Genes") {
                # a signal to use all genes instead            
                universeEvents <- NULL 
            } else if(input$universe_GO == "Selected ASE Modality") {
                if(is_valid(input$GO_EventType)) {
                    universeEvents <- res_all$EventName[
                        res_all$EventType %in% input$GO_EventType
                    ]                
                }
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

                settings_GO$resGO <- .format_GO_result(
                    .ora_internal(
                        ref(get_se())[["ontology"]], 
                        geneIds$genes, geneIds$universe,
                        ontologyType, pAdjustMethod = "BH"
                    )
                )

                settings_GO$ggplot <- .generate_ggplot_GO(settings_GO$resGO)
                settings_GO$final_plot <- ggplotly(
                    settings_GO$ggplot, tooltip = "text"
                )
            })
            settings_GO$errorMsg <- ""
        })

        output$plot_GO <- renderPlotly({
            validate(need(nrow(settings_GO$resGO) > 0, "Zero results"))
            validate(need(!is_valid(settings_GO$errorMsg), settings_GO$errorMsg))
            settings_GO$final_plot
        })

        get_ggplot <- reactive({
            settings_GO$ggplot
        })
        spModule <- vis_ggplot_server("GOplotSave", get_ggplot, volumes)

        return(settings_GO)
    })
    
}

