server_cov2 <- function(
        id, refresh_tab, volumes, get_se, get_de, get_go,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {

        settings_Cov <- setreactive_Cov2()
        
        # Reactives to NxtSE derivatives
        
        ref_r <- reactive({
            se <- get_se()
            req(is(se, "NxtSE"))
            ref(se)
        })

        colData_r <- reactive({
            se <- get_se()
            req(is(se, "NxtSE"))
            as.data.frame(colData(se))
        })
        
        rowData_r <- reactive({
            se <- get_se()
            req(is(se, "NxtSE"))
            as.data.frame(rowData(se))
        })

        # Load chromosome and gene list as soon as NxtSE is loaded
        observeEvent(ref_r(), {
            geneList <- ref_r()[["geneList"]]
            seqInfo <- ref_r()[["seqInfo"]]

            req(names(seqInfo))
            updateSelectInput(session = session, inputId = "chr_cov", 
                choices = c("(none)", names(seqInfo)), selected = "(none)")
            
            req(nrow(geneList) > 0)
            message("Populating drop-down box with ", nrow(geneList), " genes")
            updateSelectizeInput(
                session = session, inputId = "genes_cov", server = TRUE, 
                choices = c("(none)", geneList$gene_display_name), 
                selected = "(none)"
            )
        })

        observeEvent(rowData_r(), {
            rowData <- rowData_r()
            settings_Cov$event.ranges <- as.data.table(
                coord2GR(rowData$EventRegion))
            settings_Cov$event.ranges$EventName <- rowData$EventName
        })

        observeEvent(colData_r(), {
            colData <- colData_r()
            req(colData)
            condOptions <- colnames(colData)
            updateSelectInput(session = session, inputId = "condition_cov", 
                choices = c("(Individual Samples)", condOptions), 
                selected = "(Individual Samples)")
        })

        output$warning_cov <- renderText({
            validate(need(is(get_se(), "NxtSE"), 
                "Please load NxtSE first"
            ))
            
            validate(need(
                all(c("limit_start", "limit_end") %in% names(
                    settings_Cov$dataObj@args
                )), 
                "Coverage data not initialized"
            ))
            
            normEvent <- eventNorm_r()
            condName <- input$condition_cov
            condOptions <- colnames(colData_r())
            if(!(condName %in% condOptions)) condName <- NULL
            validate(need(!(!is_valid(normEvent) & is_valid(condName)), 
                "Normalization event must be selected for group plots"
            ))

            validate(need(
                all(c("limit_start", "limit_end") %in% names(
                    settings_Cov$plotObj@args
                )), 
                "Coverage plot object not yet generated"
            ))
        })
        
################### - Copied over from old dash_cov_server - ##################

        # Reactive to generate filtered DE object
        observe({
            req(get_de())

            if(input$modeFilter_COV == "Highlighted (selected) events") {
                tmpres <- as.data.table(
                    .get_unified_volcano_data(get_de()[rows_selected(),]))
            } else {
                tmpres <- as.data.table(
                    .get_unified_volcano_data(get_de()[rows_all(),]))
                    
                if(input$modeFilter_COV == "Top Gene Ontology Categories") {
                    req(get_go())
                    req(input$GOterm_COV)
                    
                    goInfo <- get_go()
                    selGOterm <- isolate(input$GOterm_COV)
                    
                    go_id <- goInfo$go_id[match(selGOterm, goInfo$Term)]
                    events <- subset_EventNames_by_GO(tmpres$EventName, go_id,
                        isolate(get_se()))
                    
                    tmpres <- tmpres[get("EventName") %in% events]
                }
            }
            
            tmpres <- tmpres[, c("EventName", "EventType"), with = FALSE]
            if(nrow(tmpres) > input$slider_num_events_cov) {
                tmpres <- tmpres[seq_len(input$slider_num_events_cov)]
            }
            
            settings_Cov$useDE <- tmpres
        })
        
        # Reactive to Populate events
        observeEvent(settings_Cov$useDE, {
            req(settings_Cov$useDE)
            res <- isolate(settings_Cov$useDE)
            
            .server_cov2_update_events_list(session, res$EventName,
                isolate(input$events_cov))
        })

        # Reactive to generate GO conditional ddb with GO terms
        observeEvent(get_go(), {
            # Update GO terms (if GO is available)
            req(get_go())
            goTerms <- isolate(get_go()$Term)
            selGOterm <- isolate(input$GOterm_COV)
            
            if(length(goTerms) > 50) goTerms <- goTerms[seq_len(50)]
            if(
                is_valid(selGOterm) && 
                selGOterm %in% goTerms
            ) {
                updateSelectInput(
                    session = session, inputId = "GOterm_COV", 
                    choices = goTerms, 
                    selected = selGOterm
                )
            } else {
                updateSelectInput(
                    session = session, inputId = "GOterm_COV", 
                    choices = goTerms
                )
            }
        })

        observeEvent(input$track_table,{
            req(input$track_table)
            settings_Cov$trackTable <- hot_to_r(input$track_table) 
        })

################### ####################################### ##################

# Triggers

# Text-entry location triggers
        chr_r <- reactive({
            req(is_valid(input$chr_cov))
            req(input$chr_cov %in% names(ref_r()$seqInfo))
            input$chr_cov
        })
        start_r <- reactive({
            req(input$start_cov)
            input$start_cov
        })
        end_r <- reactive({
            req(input$end_cov)
            input$end_cov
        })
        chr_rd <- chr_r
        start_rd <- start_r %>% debounce(1000)
        end_rd <- end_r %>% debounce(1000)

        typedGR <- eventReactive(list(
                chr_rd(), start_rd(), end_rd()
        ), {
            req(input$chr_cov)
            req(input$chr_cov %in% names(ref_r()$seqInfo))
            
            seqInfo <- ref_r()$seqInfo[input$chr_cov]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)
            
            iR <- .server_cov2_getNewIRanges(
                input$start_cov, input$end_cov, seqmax)
            
            if(length(iR) == 0) return(GRanges())
            return(GRanges(input$chr_cov, iR))
        })

# Tracks selection trigger
        tracks_r <- reactive({
            req(is_valid(settings_Cov$trackTable))

            tbl <- isolate(settings_Cov$trackTable)
            req(all(c("sample", "id") %in% colnames(tbl)))

            trackNum <- 1
            trackList <- list()            
            while(TRUE) {
                samples <- tbl$sample[which(
                    as.character(tbl$id) == as.character(trackNum)
                )]
                if(length(samples) == 0) break
                
                trackList[[trackNum]] <- as.character(samples)
                trackNum <- trackNum + 1
            }
            
            return(trackList)
        })
        tracks_rd <- tracks_r %>% debounce(3000)
        
        diff_r <- reactive({
            if(!is_valid(input$diffA) | !is_valid(input$diffB)) {
                return(list())
            }
            return(list(c(input$diffA, input$diffB)))
        })
        diff_rd <- diff_r %>% debounce(1000)

        eventNorm_r <- reactive({
            return(input$event_norm_cov)
        })
        eventNorm_rd <- eventNorm_r %>% debounce(1000)
        
        strand_r <- reactive({
            req(input$strand_cov)
            input$strand_cov
        })
        strand_rd <- strand_r %>% debounce(1000)

# Aggregate non-range trigger
        trigger <- eventReactive(list(
            input$plot_ribbon, input$plot_Jn_cov, input$normalizeCov,
            input$plot_key_iso, input$condense_cov,
            tracks_rd(), diff_rd(), eventNorm_rd(), strand_rd()
        ), {
            runif(1)
        })

# Locate by Gene or Event triggers change in chr/start/end
        eventGR <- eventReactive(input$events_cov, {
            req(input$events_cov)
            req(input$events_cov != "(none)")
            
            events_id_view <- settings_Cov$event.ranges[
                get("EventName") == input$events_cov]
            
            span        <- events_id_view$end[1] - events_id_view$start[1]
            view_start  <- max(1, events_id_view$start[1] - span)
            view_end    <- view_start + 3 * span
            
            GRanges(events_id_view$seqnames[1], IRanges(
                view_start, view_end
            ))
        })

        genesGR <- eventReactive(input$genes_cov, {
            req(input$genes_cov)
            req(input$genes_cov != "(none)")
            
            gene_id_view <- ref_r()$geneList[
                get("gene_display_name") == input$genes_cov]

            GRanges(gene_id_view$seqnames[1], IRanges(
                gene_id_view$start[1], gene_id_view$end[1]
            ))
        })

        zoomOutGR <- eventReactive(input$zoom_out_cov, {
            req(input$zoom_out_cov)
            req(input$chr_cov) 
            req(input$chr_cov %in% names(ref_r()$seqInfo))

            view_start  <- input$start_cov
            view_end    <- input$end_cov
            req(view_start, view_end, view_end - view_start >= 50)

            seqInfo <- ref_r()$seqInfo[input$chr_cov]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)

            center      <- round((view_start + view_end) / 2)
            span        <- view_end - view_start
            # zoom range is 50 * 3^z
            cur_zoom    <- floor(log(span/50) / log(3))
            new_span <- round(span * 3)
            new_start <- max(1, center - round(new_span / 2))

            GRanges(input$chr_cov, IRanges(
                new_start, new_start + new_span
            ))
        })

        zoomInGR <- eventReactive(input$zoom_in_cov, {
            req(input$zoom_in_cov)

            view_start  <- input$start_cov
            view_end    <- input$end_cov
            req(view_start, view_end, view_end - view_start >= 50)
            
            # get center of current range
            center      <- round((view_start + view_end) / 2)
            span        <- view_end - view_start
            # zoom range is 50 * 3^z
            cur_zoom    <- floor(log(span/50) / log(3))
            new_span <- round(span / 3)
            if(new_span < 50) new_span <- 50
            new_zoom <- floor(log(new_span / 50) / log(3))
            new_start <- max(1, center - round(new_span / 2))
            
            GRanges(input$chr_cov, IRanges(
                new_start, new_start + new_span
            ))
        })

        plotly_relayout <- reactive({
            req(settings_Cov$plotCount > 0)
            sourceName <- paste0("plotly_ViewRef_",
                as.character(settings_Cov$plotCount))
            event_data("plotly_relayout", source = sourceName)
        })

        plotUpdateGR <- eventReactive(plotly_relayout(), {
            layoutData <- isolate(plotly_relayout())
            chrName <- isolate(input$chr_cov)
            
            # print(layoutData)
            req(length(layoutData) == 2)
            req(all(c("xaxis.range[0]", "xaxis.range[1]") %in% 
                names(layoutData)))
            
            new_start <- max(1, round(layoutData[["xaxis.range[0]"]]))
            new_end <- round(layoutData[["xaxis.range[1]"]])
            
            # Enforce chromosome boundary
            seqInfo <- ref_r()$seqInfo[chrName]
            seqmax  <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            if(new_end > seqmax) {
                new_end <- seqmax
                if(new_end - new_start < 50) {
                    new_start <- new_end - 50
                }
            }
            # Enforce min width > 50
            if(new_end - new_start < 50) {
                if(new_end >= 50) {
                    new_start <- new_end - 50
                } else {
                    new_end <- new_start + 50
                }
            }

            GRanges(chrName, IRanges(
                new_start, new_end
            ))
        })

############################### Aggregate GRanges refresher ####################

        # function to change newGR
        updateGR <- function(gr) {
            settings_Cov$newGR <- gr
            invisible()
        }

        observeEvent(eventGR(), {
            # message("Updating event GRanges")
            updateGR(eventGR())
        })
        observeEvent(genesGR(), {
            # message("Updating gene GRanges")
            updateGR(genesGR())
        })
        observeEvent(zoomInGR(), {
            # message("Updating zoom-in GRanges")
            updateGR(zoomInGR())
        })
        observeEvent(zoomOutGR(), {
            # message("Updating zoom-out GRanges")
            updateGR(zoomOutGR())
        })
        observeEvent(plotUpdateGR(), {
            # message("Updating plot event GRanges")
            updateGR(plotUpdateGR())
        })
        observeEvent(typedGR(), {
            # message("Updating typed-in GRanges")
            updateGR(typedGR())
        })

############################### Aggregate Triggers ####################
        
        observeEvent(settings_Cov$newGR, {
            req(length(settings_Cov$newGR) == 1)
            gr <- isolate(settings_Cov$newGR)

            chrList <- names(ref_r()$seqInfo)
            updateSelectInput(session = session, inputId = "chr_cov", 
                choices = c("(none)", chrList),
                selected = as.character(seqnames(gr)))
            
            updateTextInput(session = session, inputId = "start_cov", 
                value = start(gr))
            updateTextInput(session = session, inputId = "end_cov", 
                value = end(gr))
                
            settings_Cov$plotTrigger <- runif(1)
        })
        
        observeEvent(trigger(), {
            settings_Cov$plotTrigger <- runif(1)
        })
######################### Tracks management ####################################
        
        output$track_table <- renderRHandsontable({
            .server_expr_gen_HOT(settings_Cov$trackTable, enable_select = TRUE)
        })
        
# change in condition or colData will update tracks table
        observe({
            req(get_se())
            req(input$condition_cov)
            colData <- colData_r()
            condSelected <- input$condition_cov

            trackOptions <- c()
            condOptions <- c()

            if(condSelected == "(Individual Samples)") {
                trackOptions <- rownames(colData)
            } else if(condSelected %in% colnames(colData)) {
                trackOptions <- unique(as.character(unname(unlist(
                    colData[, condSelected]))))
                condOptions <- trackOptions
            }

            updateSelectInput(session = session, inputId = "diffA", 
                choices = c("(none)", condOptions)
            )
            updateSelectInput(session = session, inputId = "diffB", 
                choices = c("(none)", condOptions)
            )
            
            df <- data.frame(
                sample = trackOptions,
                id = "", stringsAsFactors = FALSE
            )
            if(nrow(df) > 1) {
                df$id[1:2] <- c("1","2")
            } else if(nrow(df) == 1) {
                df$id <- "1"
            }
            settings_Cov$trackTable <- df
        })

############################  trigger ####################################

        observeEvent(settings_Cov$plotTrigger, {
            req(get_se())
            req(all(is_valid(covfile(get_se()))))

            gr <- isolate(settings_Cov$newGR)
            req(length(gr) > 0)

            condName <- isolate(input$condition_cov)
            req(condName)

            trackList <- isolate(tracks_r())
            diffList <- isolate(diff_r())
            req(trackList)

        # Grab everything we need from the start (except normEvent)
            tmpChr <- as.character(seqnames(gr))
            tmpStart <- start(gr)
            tmpEnd <- end(gr)

            colData <- isolate(colData_r())
            condOptions <- colnames(colData)
            if(length(condOptions) == 0) {
                condName <- NULL
            } else if(!(condName %in% condOptions)) {
                condName <- NULL
            }

            trackOptions <- NULL
            if(is.null(condName)) {
                condName <- NULL
                trackOptions <- rownames(colData)
            } else {
                trackOptions <- unique(as.character(
                    unname(unlist(colData[,condName]))))
            }
            strand <- isolate(input$strand_cov)

        # Starting deal-breakers
            req(length(trackOptions) > 0)

        # Do we need to update cDO? If out of range, or if no cDO
            refreshCDO <- FALSE
            args <- settings_Cov$dataObj@args
            if(!all(c("limit_start", "limit_end") %in% names(args))) {
                # Likely not a valid cDO, regenerate it
                refreshCDO <- TRUE
            } else if(
                    args[["view_chr"]] != tmpChr ||
                (
                    args[["limit_start"]] > tmpStart |
                    args[["limit_end"]] < tmpEnd
                )
            ) {
                refreshCDO <- TRUE  
            }

            dataObj <- isolate(settings_Cov$dataObj)
            if(refreshCDO) {
                withProgress(message = 'Retrieving COV data...', value = 0, {
                    dataObj <- getCoverageData(
                        isolate(get_se()),
                        seqname = tmpChr,
                        start = tmpStart,
                        end = tmpEnd,
                        tracks = colnames(isolate(get_se()))
                    )
                })
            }

        # Check cDO is valid before assigning to reactive list
            args <- isolate(dataObj@args)
            req(all(c("limit_start", "limit_end") %in% names(args)))
            settings_Cov$dataObj <- dataObj

        # Update Norm Event options and retrieve normEvent
            normEvent <- isolate(eventNorm_r())
            normRowData <- isolate(dataObj@normData$rowData)
            availEvents <- c()
            if(nrow(normRowData) > 0) {
                availEvents <- normRowData$EventName
                requestedEvent <- isolate(input$events_cov)
                if(
                        requestedEvent %in% availEvents &
                        !(normEvent %in% availEvents)
                ) {
                    normEvent <- requestedEvent
                } else if(!(normEvent %in% availEvents)) {
                    normEvent <- "(none)"
                }
            } else {
                normEvent <- "(none)"
            }
            updateSelectInput(session = session, inputId = "event_norm_cov", 
                choices = c("(none)", availEvents),
                selected = normEvent
            )

            # Abort if valid condName but no valid normEvent
            req(!(
                !is_valid(normEvent) & is_valid(condName)
            ))
            
            # Do we need to update cPO?
            # - change in condition (which affects tracks), 
            # - strand, or normalization event
 
            refreshCPO <- refreshCDO            
            args <- isolate(settings_Cov$plotObj@args)
            if(!("tracks" %in% names(args))) {
                # Invalid cPO
                refreshCPO <- TRUE  
            } else if(args[["strand"]] != strand) {
                # Different strand requested
                refreshCPO <- TRUE            
            } else if(!all(trackOptions %in% args[["tracks"]])) {
                # Some tracks missing from cPO
                refreshCPO <- TRUE  
            } else if((!is_valid(normEvent) && "Event" %in% names(args))) {
                # Need to remove `Event` from current cPO
                refreshCPO <- TRUE              
            } else if((is_valid(normEvent) && !("Event" %in% names(args)))) {
                # Need to assign `Event` to a cPO that doesn't have one
                refreshCPO <- TRUE              
            } else if((is_valid(normEvent) && args[["Event"]] != normEvent)) {
                # Need to assign a different event to current cPO
                refreshCPO <- TRUE              
            }

            plotObj <- isolate(settings_Cov$plotObj)
            if(refreshCPO) {
                withProgress(
                    message = 'Calculating track coverages...', 
                    value = 0, 
                {
                    cPOargList <- list(
                        object = dataObj,
                        strand = strand,
                        tracks = trackOptions
                    )
                    if(is_valid(condName)) cPOargList[["condition"]] <- condName
                    if(is_valid(normEvent)) cPOargList[["Event"]] <- normEvent
                    plotObj <- do.call(getPlotObject, cPOargList)
                })
            }
            
            # Check if cPO is valid before assigning to reactive list
            args <- plotObj@args
            req(all(c("limit_start", "limit_end") %in% names(args)))
            settings_Cov$plotObj <- plotObj

            refreshPlotly <- refreshCPO
            oldPlotSettings <- isolate(settings_Cov$oldPlotSettings)
            newSettings <- list(
                view_start = tmpStart, view_end = tmpEnd,
                trackList = trackList, diffList = diffList,
                ribbon_mode = isolate(input$plot_ribbon),
                plotJunctions = isolate(input$plot_Jn_cov),
                normalizeCoverage = isolate(input$normalizeCov),
                filterByEventTranscripts = isolate(input$plot_key_iso),
                condenseTranscripts = isolate(input$condense_cov)
            )
            if(!identical(newSettings, oldPlotSettings)) {
                refreshPlotly <- TRUE
            }
            # TODO - plot annotation track only
            if(length(trackList) == 0) refreshPlotly <- FALSE

            if(refreshPlotly) {
                withProgress(
                    message = 'Generating plot...', 
                    value = 0, 
                {
                    plotlyObj <- plotView(
                        plotObj, oldP = isolate(settings_Cov$plotlyObj),
                        view_start = newSettings[["view_start"]], 
                        view_end = newSettings[["view_end"]],
                        trackList = newSettings[["trackList"]],
                        diffList = newSettings[["diffList"]],
                        diff_stat = "t-test",
                        ribbon_mode = newSettings[["ribbon_mode"]],
                        plotJunctions = newSettings[["plotJunctions"]],
                        normalizeCoverage = newSettings[["normalizeCoverage"]],
                        filterByEventTranscripts = newSettings[["filterByEventTranscripts"]],
                        condenseTranscripts = newSettings[["condenseTranscripts"]],
                        usePlotly = TRUE
                    )
                })
            } else {
                plotlyObj <- NULL
            }
            
            req(is(plotlyObj, "covPlotly"))
            settings_Cov$oldPlotSettings <- newSettings
            settings_Cov$exons_gr <- getExonRanges(plotlyObj)            
            settings_Cov$plotlyObj <- setResolution(plotlyObj, 
                isolate(input$slider_num_plotRes))
                
            # Increment plot count to trigger synthFig()
            settings_Cov$plotCount <- isolate(settings_Cov$plotCount) + 1
        })

        # Trigger new plotlyFig every time this increments
        synthFig <- eventReactive(settings_Cov$plotCount, {
            req(settings_Cov$plotlyObj)
            req(all(
                c("xrange", "resolution") %in% 
                names(settings_Cov$plotlyObj@args)
            ))
            plotCount <- isolate(settings_Cov$plotCount)

            doExons <- isolate(input$exonMode_cov)
            if(doExons) {
                fig <- .covPlotlyMake(settings_Cov$plotlyObj, showExons = TRUE)            
            } else {
                fig <- .covPlotlyMake(settings_Cov$plotlyObj)
            }

            fig$x$source <- paste0("plotly_ViewRef_",
                as.character(plotCount))
            if(packageVersion("plotly") >= "4.9.0") {
                event_register(fig, "plotly_relayout")
            }
            
            return(fig)
        })

        # If this changes and plotlyObj is a valid covPlotly, increment by 1
        observeEvent(input$exonMode_cov, {
            req(settings_Cov$plotlyObj)
            req(all(
                c("xrange", "resolution") %in% 
                names(settings_Cov$plotlyObj@args)
            ))
            settings_Cov$plotCount <- isolate(settings_Cov$plotCount) + 1
        })

        observeEvent(input$slider_num_plotRes, {
            req(settings_Cov$plotlyObj)
            req(all(
                c("xrange", "resolution") %in% 
                names(settings_Cov$plotlyObj@args)
            ))
            req(
                settings_Cov$plotlyObj@args[["resolution"]] !=
                input$slider_num_plotRes
            )
            settings_Cov$plotlyObj <- setResolution(
                settings_Cov$plotlyObj, input$slider_num_plotRes
            )
            
            settings_Cov$plotCount <- isolate(settings_Cov$plotCount) + 1
        })
        
        output$plot_cov <- renderPlotly({
            synthFig()
        })

        observeEvent(settings_Cov$exons_gr, {
            gr <- isolate(settings_Cov$exons_gr)
            if(length(gr) > 0) {
                df <- data.frame(
                    exon = sort(names(gr)),
                    selected = FALSE
                )
            } else {
                df <- data.frame()
            }
            settings_Cov$exonsTable <- df
        })

        output$exons_lookup <- renderRHandsontable({
            .server_expr_gen_HOT(
                settings_Cov$exonsTable, 
                enable_select = TRUE,
                lockedColumns = "exon"
            )
        })

        observeEvent(input$exons_lookup,{
            req(input$exons_lookup)
            settings_Cov$exonsTable <- hot_to_r(input$exons_lookup) 
        })

        exonsSelected_r <- reactive({
            df_exons <- settings_Cov$exonsTable
            req(df_exons)
            req(nrow(df_exons) > 0)
            
            if(sum(df_exons$selected) < 2) return(NULL)
            return(df_exons$exon[df_exons$selected == TRUE])
        })
        exonsSelected_rd <- exonsSelected_r %>% debounce(1000)

        synthExonsFig <- eventReactive(exonsSelected_rd(), {
            req(settings_Cov$plotlyObj)
            req(all(
                c("xrange", "resolution") %in% 
                names(settings_Cov$plotlyObj@args)
            ))
        
            exonsNames <- exonsSelected_r()
        
            validate(need(exonsNames,
                "Select two or more exons to plot exon-centric coverage plot"
            ))

            gr <- isolate(settings_Cov$exons_gr)
            newSettings <- isolate(settings_Cov$oldPlotSettings)
            # print(gr[exonsNames])
            ggp <- plotView(
                isolate(settings_Cov$plotObj), 
                # view_start = newSettings[["view_start"]], 
                # view_end = newSettings[["view_end"]],
                trackList = newSettings[["trackList"]],
                diffList = newSettings[["diffList"]],
                diff_stat = "t-test",
                ribbon_mode = newSettings[["ribbon_mode"]],
                plotJunctions = newSettings[["plotJunctions"]],
                normalizeCoverage = newSettings[["normalizeCoverage"]],
                filterByEventTranscripts = newSettings[["filterByEventTranscripts"]],
                condenseTranscripts = newSettings[["condenseTranscripts"]],
                plotRanges = gr[exonsNames],
                usePlotly = FALSE
            )

            return(ggp)
        })

        output$stillplot_cov <- renderPlot({
            synthExonsFig()
        })

        return(settings_Cov)
    })
}

# Populate drop-down for searching by Event
.server_cov2_update_events_list <- function(
    session, EventNames, selectedEvent
) {
    if(selectedEvent %in% EventNames) {
        updateSelectizeInput(session = session, 
            inputId = "events_cov", server = TRUE,
            choices = c("(none)", EventNames), 
            selected = selectedEvent
        )
    } else {
        updateSelectizeInput(session = session, 
            inputId = "events_cov", server = TRUE,
            choices = c("(none)", EventNames),
            selected = "(none)"
        )
    }
}

# Error-checking for:
# - typed-in start / end coordinates
# - changing chromosome at drop-down box
.server_cov2_getNewIRanges <- function(target_start, target_end, seqmax) {
    # do nothing if empty 
    if(!is_valid(target_start) | !is_valid(target_end)) return(IRanges())
    
    # do nothing if not numeric
    if(is_valid(target_start) && !is.numeric(target_start))  return(IRanges())
    if(is_valid(target_end) && !is.numeric(target_end))  return(IRanges())
    
    # Cap target_end at seqmax
    if(target_end > seqmax) target_end <- seqmax
    
    # Cap min width at 50 - this assumes seqlength > 50
    if(target_end < 50) target_end <- 50
    if(target_end - target_start < 50) target_start <- target_end - 50

    return(IRanges(target_start, target_end))
}
