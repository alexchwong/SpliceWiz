server_cov2 <- function(
        id, refresh_tab, volumes, get_se, get_de, get_go,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {

        settings_Cov <- setreactive_Cov2()
        get_ref <- reactive({
            se <- get_se()
            req(is(se, "NxtSE"))
            ref(se)
        })
        
        observeEvent(refresh_tab(), {
            req(refresh_tab())
            output$warning_cov <- renderText({
                validate(need(is(get_se(), "NxtSE"), 
                    "Please build experiment first"))
            })
            
            se <- isolate(get_se())
            req(is(se, "NxtSE"))
            
            # Slow step - populate genes
            ref <- isolate(get_ref())
            req(ref)
            
            withProgress(message = 'Populating gene list...', value = 0, {
                .server_cov2_update_genes(session, ref$geneList)
            })
            
            settings_Cov$event.ranges <- as.data.table(
                coord2GR(rowData(se)$EventRegion))
            settings_Cov$event.ranges$EventName <- rowData(se)$EventName
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
            req(input$chr_cov %in% names(get_ref()$seqInfo))
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
            req(input$chr_cov %in% names(get_ref()$seqInfo))
            
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
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
            
            gene_id_view <- get_ref()$geneList[
                get("gene_display_name") == input$genes_cov]

            GRanges(gene_id_view$seqnames[1], IRanges(
                gene_id_view$start[1], gene_id_view$end[1]
            ))
        })

        zoomOutGR <- eventReactive(input$zoom_out_cov, {
            req(input$zoom_out_cov)
            req(input$chr_cov) 
            req(input$chr_cov %in% names(get_ref()$seqInfo))

            view_start  <- input$start_cov
            view_end    <- input$end_cov
            req(view_start, view_end, view_end - view_start >= 50)


            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)

            center      <- round((view_start + view_end) / 2)
            span        <- view_end - view_start
            # zoom range is 50 * 3^z
            cur_zoom    <- floor(log(span/50) / log(3))

            new_span <- round(span * 3)
            # if(new_span > seqmax - 1) new_span = seqmax - 1
            new_start <- max(1, center - round(new_span / 2))
            
            cur_zoom <- floor(log(new_span/50) / log(3))

            GRanges(input$chr_cov, IRanges(
                new_start, new_start + new_span
            ))
        })

        zoomInGR <- eventReactive(input$zoom_out_cov, {
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

            cur_zoom <- floor(log(new_span/50) / log(3))
            # output$label_zoom_cov <- renderText({16 - cur_zoom})
            
            GRanges(input$chr_cov, IRanges(
                new_start, new_start + new_span
            ))
        })

        plotly_relayout <- reactive({
            req(length(settings_Cov$oldPlotSettings) > 0)
            req(settings_Cov$plot_ini == TRUE)
            event_data("plotly_relayout", source = "plotly_ViewRef")
        })

        plotUpdateGR <- eventReactive(plotly_relayout(), {
            layoutData <- isolate(plotly_relayout())

            print(layoutData)
            req(length(layoutData) == 2)
            req(all(c("xaxis.range[0]", "xaxis.range[1]") %in% 
                names(layoutData)))
            
            new_start <- max(1, round(layoutData[["xaxis.range[0]"]]))
            new_end <- round(layoutData[["xaxis.range[1]"]])
            
            # Enforce chromosome boundary
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
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

            GRanges(input$chr_cov, IRanges(
                new_start, new_end
            ))
        })

############################### Aggregate GRanges refresher ####################

        observeEvent(eventGR(), {
            settings_Cov$newGR <- eventGR()
        })
        observeEvent(genesGR(), {
            settings_Cov$newGR <- genesGR()
        })
        observeEvent(zoomInGR(), {
            settings_Cov$newGR <- zoomInGR()
        })
        observeEvent(zoomOutGR(), {
            settings_Cov$newGR <- zoomOutGR()
        })
        observeEvent(plotUpdateGR(), {
            settings_Cov$newGR <- plotUpdateGR()
        })
        observeEvent(typedGR(), {
            settings_Cov$newGR <- typedGR()
        })

############################### Aggregate Triggers ####################
        
        observeEvent(settings_Cov$newGR, {
            gr <- isolate(settings_Cov$newGR)
            chrList <- names(get_ref()$seqInfo)
            updateSelectInput(session = session, inputId = "chr_cov", 
                choices = c("(none)", chrList),
                selected = seqnames(gr))
            
            updateTextInput(session = session, inputId = "start_cov", 
                value = start(gr))
            updateTextInput(session = session, inputId = "end_cov", 
                value = start(gr))
                
            settings_Cov$plotTrigger <- runif(1)
        })
        
        observeEvent(trigger(), {
            settings_Cov$plotTrigger <- runif(1)
        })
######################### Tracks management ####################################
        
        output$track_table <- renderRHandsontable({
            .server_expr_gen_HOT(settings_Cov$trackTable, enable_select = TRUE)
        })
        
# colData update changes condition drop-down box
        observeEvent(get_se(), {
            colData <- isolate(colData(get_se()))

            conditions <- colnames(colData)
            updateSelectInput(session = session, inputId = "condition_cov", 
                choices = c("(Individual Samples)", colnames(colData))
            )
            df <- data.frame(
                sample = colnames(get_se()),
                id = "", stringsAsFactors = FALSE
            )
            if(nrow(df) > 1) {
                df$id[1:2] <- c("1","2")
            } else if(nrow(df) == 1) {
                df$id <- "1"
            }
            settings_Cov$trackTable <- df
        })

# change in condition will update tracks table
        observeEvent(input$condition_cov, {
            req(input$condition_cov)
            req(get_se())
            colData <- as.data.frame(isolate(colData(get_se())))

            trackOptions <- c()
            condOptions <- c()

            if(input$condition_cov == "(Individual Samples)") {
                trackOptions <- colnames(get_se())
            } else if(input$condition_cov %in% colnames(colData)) {
                trackOptions <- unique(as.character(unname(unlist(
                    colData[, input$condition_cov]))))
                condOptions <- trackOptions
            }

            updateSelectInput(session = session, inputId = "diffA", 
                choices = c("(none)", trackOptions)
            )
            updateSelectInput(session = session, inputId = "diffB", 
                choices = c("(none)", trackOptions)
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

        makePlot <- eventReactive(settings_Cov$plotTrigger, {
            req(length(settings_Cov$newGR) > 0)

            tmpChr <- as.character(seqnames(settings_Cov$newGR))
            tmpStart <- start(settings_Cov$newGR)
            tmpEnd <- end(settings_Cov$newGR)

            req(get_se())
            req(input$condition_cov)
            req(is_valid(settings_Cov$trackTable))

            # Do we need to update cDO?
            refreshCDO <- FALSE
            args <- settings_Cov$dataObj@args
            if(!all(c("limit_start", "limit_end") %in% names(args))) {
                # Likely not a valid cDO, regenerate it
                refreshCDO <- TRUE
            } else {                
                if(
                        args[["view_chr"]] != tmpChr ||
                        (
                            args[["limit_start"]] > tmpStart |
                            args[["limit_end"]] < tmpEnd
                        )
                ) {
                    refreshCDO <- TRUE  
                }
            }

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
            } else {
                dataObj <- isolate(settings_Cov$dataObj)
            }

        # Check cDO is valid
            args <- isolate(dataObj@args)
            output$warning_cov <- renderText({
                validate(need(
                    all(c("limit_start", "limit_end") %in% names(args)), 
                    "covDataObject invalid"
                ))
            })
            req(all(c("limit_start", "limit_end") %in% names(args)))
            settings_Cov$dataObj <- dataObj

        # Update Norm Event options
            normData <- isolate(settings_Cov$dataObj@normData$rowData)
            normEvent <- "(none)"
            availEvents <- c()
            if(nrow(normData) > 0) {
                availEvents <- normData$EventName
                requestedEvent <- isolate(input$events_cov)
                normEvent <- isolate(eventNorm_rd())
                
                if(
                        requestedEvent %in% availEvents &
                        !(normEvent %in% availEvents)
                ) {
                    normEvent <- requestedEvent
                } else if(!(normEvent %in% availEvents)) {
                    normEvent <- "(none)"
                }
            }
            updateSelectInput(session = session, inputId = "event_norm_cov", 
                choices = c("(none)", availEvents),
                selected = normEvent
            )

            # Do we need to update cPO?
            # - change in condition (which affects tracks), 
            # - strand, or normalization event
 
            refreshCPO <- refreshCDO            

            args <- isolate(settings_Cov$plotObj@args)
            colData <- isolate(as.data.frame(colData(get_se())))
            condName <- isolate(input$condition_cov)
            trackOptions <- NULL
            if(!(condName %in% colnames(colData))) {
                condName <- NULL
                trackOptions <- rownames(colData)
            } else {
                trackOptions <- unique(as.character(
                    unname(unlist(colData[,condName]))))
            }
            strand <- isolate(input$strand_cov)
            
            # Only update if there is a valid trackList
            if(length(trackOptions) > 0) {
                if(!("tracks" %in% names(args))) {
                    # null CPO
                    refreshCPO <- TRUE
                } else if(!all(trackOptions %in% args[["tracks"]])) {
                    refreshCPO <- TRUE
                }            
            }
            # Check normalization event
            if(!refreshCPO) {
                if(!is_valid(normEvent) && "Event" %in% names(args)) {
                    refreshCPO <- TRUE
                } else if(is_valid(normEvent) && !("Event" %in% names(args))) {
                    refreshCPO <- TRUE
                } else if(is_valid(normEvent) && args[["Event"]] != normEvent) {
                    refreshCPO <- TRUE
                }
            }
            # Check strand
            if(!refreshCPO) {
                if(args[["strand"]] != strand) {
                    refreshCPO <- TRUE
                }
            }
            
            warningTxt <- ""
            if(refreshCPO) {
                output$warning_cov <- renderText({
                    validate(need(!is_valid(normEvent) & is_valid(condName), 
                        "Norm event must be selected for group plots"))
                })
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
                    if(!is_valid(normEvent) & is_valid(condName)) {
                        plotObj <- covPlotObject()
                    } else {
                        plotObj <- do.call(getPlotObject, cPOargList)
                    }
                })
            } else {
                plotObj <- isolate(settings_Cov$plotObj)
            }

            # Check if cPO is valid
            args <- plotObj@args
            output$warning_cov <- renderText({
                validate(need(
                    all(c("limit_start", "limit_end") %in% names(args)), 
                    "covPlotObject invalid"
                ))
            })
            req(all(c("limit_start", "limit_end") %in% names(args)))
            settings_Cov$plotObj <- plotObj

            trackList <- isolate(tracks_r())
            diffList <- isolate(diff_r())

            refreshPlotly <- refreshCPO
            newSettings <- list(
                view_start = tmpStart, view_end = tmpEnd,
                trackList = trackList, diffList = diffList,
                ribbon_mode = isolate(input$plot_ribbon),
                plotJunctions = isolate(input$plot_Jn_cov),
                normalizeCoverage = isolate(input$normalizeCov),
                filterByEventTranscripts = isolate(input$plot_key_iso),
                condenseTranscripts = isolate(input$condense_cov)
            )
            if(!identical(newSettings, settings_Cov$oldPlotSettings)) {
                refreshPlotly <- TRUE
            }
            if(length(trackList) == 0) refreshPlotly <- FALSE

            if(refreshPlotly) {
                settings_Cov$plot_ini <- FALSE
                plotlyObj <- plotView(
                    plotObj, oldP = isolate(settings_Cov$plotlyObj),
                    view_start = tmpStart, view_end = tmpEnd,
                    trackList = trackList,
                    diffList = diffList,
                    diff_stat = "t-test",
                    ribbon_mode = isolate(input$plot_ribbon),
                    plotJunctions = isolate(input$plot_Jn_cov),
                    normalizeCoverage = isolate(input$normalizeCov),
                    filterByEventTranscripts = isolate(input$plot_key_iso),
                    condenseTranscripts = isolate(input$condense_cov),
                    usePlotly = TRUE
                )
                settings_Cov$oldPlotSettings <- newSettings
            } else {
                plotlyObj <- NULL
            }
            
            req(is(plotlyObj, "covPlotly"))
            settings_Cov$plotlyObj <- plotlyObj
            fig <- .covPlotlyMake(plotlyObj)
            req(fig)
            
            fig$x$source <- "plotly_ViewRef"
            if(packageVersion("plotly") >= "4.9.0") {
                event_register(fig, "plotly_relayout")
            }
            settings_Cov$plot_ini <- TRUE
            return(fig)
        })
        
        output$plot_cov <- renderPlotly({
            fig <- makePlot()
            req(fig)
            fig
        })

        return(settings_Cov)
    })
}

# Update gene drop-down
.server_cov2_update_genes <- function(
    session, geneList
) {
    if(!is.null(geneList)) {
        message("Populating drop-down box with ", 
            length(unique(geneList$gene_display_name)), " genes")
        updateSelectInput(session = session, inputId = "chr_cov", 
            choices = c("(none)", 
                as.character(sort(unique(geneList$seqnames)))),
            selected = "(none)")
        # Initialize first then repopulate later
        updateSelectizeInput(session = session, inputId = "genes_cov",
            server = TRUE, choices = c("(none)"), 
            selected = "(none)")
        updateSelectizeInput(session = session, inputId = "genes_cov",
            server = TRUE, choices = c("(none)", geneList$gene_display_name), 
            selected = "(none)")
    } else {
        updateSelectInput(session = session, inputId = "chr_cov", 
            choices = c("(none)"), selected = "(none)")
        updateSelectizeInput(session = session, inputId = "genes_cov", 
            server = TRUE, choices = c("(none)"), selected = "(none)") 
    }
}

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

# Changes start and end coordinates, if over seqmax
# also if changing to shorter chromosome and prev end > chrom length
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

.server_cov2_trackList_names2ID <- function(trackList, args) {
    if(!("tracks" %in% names(args))) {
        return(NULL)
    }
    tracks <- args[["tracks"]]
    outList <- list()
    for(i in seq_len(length(trackList))) {
        trackVec <- trackList[[i]]
        trackID <- match(trackVec, tracks)
        if(any(is.na(trackID))) return(NULL)
        outList[[i]] <- trackID
    }
    return(outList)
}

.server_cov2_trackList2Vec <- function(trackList) {
    vec <- c()
    for(i in seq_len(length(trackList))) {
        vec <- c(vec, trackList[[i]])
    }
    return(vec)
}