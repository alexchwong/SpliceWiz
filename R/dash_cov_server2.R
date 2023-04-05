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

        observeEvent(list(
            start_rd(), end_rd()
        ), {
            settings_Cov$new_range <- IRanges(start_r(), end_r())
        })

# Plot options triggers
        observeEvent(list(
            input$plot_ribbon, input$plot_Jn_cov, input$normalizeCov,
            input$plot_key_iso, input$condense_cov
        ), {
            settings_Cov$trigger <- runif(1)
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
        
        observeEvent(list(
            tracks_rd(), diff_rd(), eventNorm_rd(), strand_rd()
        ), {
            settings_Cov$trigger <- runif(1)
        })

# Start/End sanity check when these values are changed; also updates zoom
        observeEvent(chr_rd(), {
            seqInfo <- get_ref()$seqInfo[chr_rd()]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)
            output <- .server_cov2_check_start_end(
                input, session, output, seqmax)
        })
        observeEvent(list(start_rd(), end_rd()), {
            req(input$chr_cov, input$chr_cov %in% names(get_ref()$seqInfo))
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)
            output <- .server_cov2_check_start_end(
                input, session, output, seqmax)
        })

# Locate by Gene or Event triggers change in chr/start/end
        observeEvent(input$events_cov, {
            req(input$events_cov)
            req(input$events_cov != "(none)")
            
            events_id_view <- settings_Cov$event.ranges[
                get("EventName") == input$events_cov]
            .server_cov2_locate_events(input, session, events_id_view)
        })
        observeEvent(input$genes_cov, {
            req(input$genes_cov)
            req(input$genes_cov != "(none)")
            
            gene_id_view <- get_ref()$geneList[
                get("gene_display_name") == input$genes_cov]
            .server_cov2_locate_genes(input, session, gene_id_view)
        })

# Sets start and end when zooming in / out
        observeEvent(input$zoom_out_cov, {
            req(input$zoom_out_cov, input$chr_cov) 
            req(input$chr_cov %in% names(get_ref()$seqInfo))
            
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            output <- .server_cov2_zoom_out(
                input, output, session, seqInfo, settings_Cov)
        })
        observeEvent(input$zoom_in_cov, {
            req(input$zoom_in_cov)
            
            output <- .server_cov2_zoom_in(input, output, session, settings_Cov)
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

            if(input$condition_cov == "(Individual Samples)") {
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
            } else if(input$condition_cov %in% colnames(colData)) {
                condOptions <- unique(as.character(unname(unlist(
                    colData[, input$condition_cov]))))
                df <- data.frame(
                    sample = condOptions,
                    id = "", stringsAsFactors = FALSE
                )
                if(nrow(df) > 1) {
                    df$id[1:2] <- c("1","2")
                } else if(nrow(df) == 1) {
                    df$id <- "1"
                }
                settings_Cov$trackTable <- df
                updateSelectInput(session = session, inputId = "diffA", 
                    choices = c("(none)", condOptions)
                )
                updateSelectInput(session = session, inputId = "diffB", 
                    choices = c("(none)", condOptions)
                )
            }
        })

############################  trigger ####################################

        observeEvent(list(
                settings_Cov$trigger, 
                settings_Cov$new_range
        ), {
            # message("New cDO/CPO event triggered")
            req(length(settings_Cov$new_range) > 0)
            req(start(settings_Cov$new_range) > 0)
            req(end(settings_Cov$new_range) > start(settings_Cov$new_range))

            tmpStart <- start(settings_Cov$new_range)
            tmpEnd <- end(settings_Cov$new_range)
            updateTextInput(session = session, inputId = "start_cov", 
                value = tmpStart)
            updateTextInput(session = session, inputId = "end_cov", 
                value = tmpEnd)

            # - Locale-based trigger
            req(input$chr_cov)
            # Check validity of chr/start/end call
            tmpChr <- isolate(input$chr_cov)

            req(isolate(get_se()))
            req(isolate(input$condition_cov))
            req(is_valid(isolate(settings_Cov$trackTable)))

            args <- isolate(settings_Cov$dataObj@args)
            
            # Do we need to update cDO?
            refreshCDO <- FALSE
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
                # message("cDO retrieved")
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
            if(nrow(normData) == 0) {
                updateSelectInput(session = session, inputId = "event_norm_cov", 
                    choices = c("(none)"))
                normEvent <- "(none)"
            } else {
                availEvents <- normData$EventName
                requestedEvent <- isolate(input$events_cov)
                normEvent <- isolate(eventNorm_rd())
                
                if(
                        requestedEvent %in% availEvents &
                        !(normEvent %in% availEvents)
                ) {
                    updateSelectInput(
                        session = session, inputId = "event_norm_cov", 
                        choices = c("(none)", normData$EventName),
                        selected = requestedEvent
                    )
                    normEvent <- requestedEvent
                } else if(normEvent %in% availEvents) {
                    updateSelectInput(
                        session = session, inputId = "event_norm_cov", 
                        choices = c("(none)", normData$EventName),
                        selected = normEvent
                    )
                } else {
                    updateSelectInput(
                        session = session, inputId = "event_norm_cov", 
                        choices = c("(none)", normData$EventName)
                    )
                    normEvent <- "(none)"
                }
            }

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
            event <- NULL
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
                    if(is_valid(normEvent) & is_valid(condName)) {
                        plotObj <- getPlotObject(
                            dataObj,
                            Event = normEvent,
                            strand = strand,
                            tracks = trackOptions,
                            condition = condName
                        )
                        # message("cPO retrieved")
                    } else if(is_valid(normEvent) & !is_valid(condName)) {
                        plotObj <- getPlotObject(
                            dataObj, 
                            Event = normEvent,
                            strand = strand,
                            tracks = trackOptions
                        )
                        # message("cPO retrieved")
                    } else if(!is_valid(normEvent) & is_valid(condName)) {
                        # condition is selected before event name
                        # fire a warning and fail to get valid covPlotObject
                        plotObj <- covPlotObject()
                    } else {
                        plotObj <- getPlotObject(
                            dataObj,
                            strand = strand,
                            tracks = trackOptions
                        )
                        # message("cPO retrieved")
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
                # message("covPlotly retrieved")
            } else {
                plotlyObj <- NULL
            }

            req(is(plotlyObj, "covPlotly"))
            settings_Cov$plotlyObj <- plotlyObj
            fig <- .covPlotlyMake(plotlyObj)
            req(fig)
            
            output$plot_cov <- renderPlotly({
                req(fig)
                
                fig$x$source <- "plotly_ViewRef"
                if(packageVersion("plotly") >= "4.9.0") {
                    event_register(fig, "plotly_relayout")
                }
                print(fig)
            })
        })
        
        # Allow update locale on zoom / pan
        settings_Cov$plotly_relayout <- reactive({
            req(length(newSettings) > 0)
            event_data("plotly_relayout", source = "plotly_ViewRef")
        })
        observeEvent(settings_Cov$plotly_relayout(), {
            layoutData <- isolate(settings_Cov$plotly_relayout())
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
          
            updateTextInput(session = session, inputId = "start_cov", 
                value = new_start)
            updateTextInput(session = session, inputId = "end_cov", 
                value = new_end)
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
.server_cov2_check_start_end <- function(input, session, output, seqmax) {
    target_start    <- isolate(input$start_cov)
    target_end      <- isolate(input$end_cov)
    # do nothing if empty 
    if(!is_valid(target_start) | !is_valid(target_end)) return(output)
    
    # do nothing if not numeric
    if(is_valid(target_start) && !is.numeric(target_start))  return(output)
    if(is_valid(target_end) && !is.numeric(target_end))  return(output)
    
    # Cap target_end at seqmax
    if(target_end > seqmax) target_end <- seqmax
    
    # Cap min width at 50 - this assumes seqlength > 50
    if(target_end < 50) target_end <- 50
    if(target_end - target_start < 50) target_start <- target_end - 50
    
    # Work out the zoom-factor
    span <- target_end - target_start
    cur_zoom <- floor(log(span/50) / log(3))
    
    output$label_zoom_cov <- renderText({16 - cur_zoom})
    if(target_end != input$end_cov)
        updateTextInput(session = session, inputId = "end_cov", 
            value = target_end)
    if(target_start != input$start_cov)
    updateTextInput(session = session, inputId = "start_cov", 
        value = target_start)
    return(output)
}

# Sets start and end given an Event
.server_cov2_locate_events <- function(input, session, events_id_view) {
        
    span        <- events_id_view$end[1] - events_id_view$start[1]
    view_start  <- max(1, events_id_view$start[1] - span)
    view_end    <- view_start + 3 * span
    
    updateSelectInput(session = session, inputId = "chr_cov", 
        selected = events_id_view$seqnames[1])
    updateTextInput(session = session, inputId = "start_cov", 
        value = view_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = view_end)
}

# Sets start and end given a Gene
.server_cov2_locate_genes <- function(input, session, gene_id_view) {
    updateSelectInput(session = session, inputId = "chr_cov", 
        selected = gene_id_view$seqnames[1])
    updateTextInput(session = session, inputId = "start_cov", 
        value = gene_id_view$start[1])
    updateTextInput(session = session, inputId = "end_cov", 
        value = gene_id_view$end[1])
}

# Zoom out
.server_cov2_zoom_out <- function(
        input, output, session, seqInfo, settings_Cov
) {
    view_start  <- input$start_cov
    view_end    <- input$end_cov
    req(view_start, view_end, view_end - view_start >= 50)
    
    seqmax      <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
    # get center of current range
    center      <- round((view_start + view_end) / 2)
    span        <- view_end - view_start
    # zoom range is 50 * 3^z
    cur_zoom    <- floor(log(span/50) / log(3))

    new_span <- round(span * 3)
    # if(new_span > seqmax - 1) new_span = seqmax - 1
    new_start <- max(1, center - round(new_span / 2))
    
    settings_Cov$plot_params$start <- new_start
    settings_Cov$plot_params$end <- new_start + new_span
    cur_zoom <- floor(log(new_span/50) / log(3))
    output$label_zoom_cov <- renderText({16 - cur_zoom})
    
    updateTextInput(session = session, inputId = "start_cov", 
        value = new_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = new_start + new_span)
    return(output)
}

# Zoom in
.server_cov2_zoom_in <- function(input, output, session, settings_Cov) {
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

    settings_Cov$plot_params$start <- new_start
    settings_Cov$plot_params$end <- new_start + new_span
    cur_zoom <- floor(log(new_span/50) / log(3))
    output$label_zoom_cov <- renderText({16 - cur_zoom})
    
    updateTextInput(session = session, inputId = "start_cov", 
        value = new_start)
    updateTextInput(session = session, inputId = "end_cov", 
        value = new_start + new_span)
    return(output)
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