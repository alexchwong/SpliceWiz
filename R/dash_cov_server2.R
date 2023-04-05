server2_cov <- function(
        id, refresh_tab, volumes, get_se, get_de, get_go,
        rows_all, rows_selected
) {
    moduleServer(id, function(input, output, session) {

        settings_Cov <- setreactive_Cov2()
        get_ref <- reactive{
            se <- get_se()
            ref(se)
        }
        
        observeEvent(refresh_tab(), {
            req(refresh_tab())
            output$warning_cov <- renderText({
                validate(need(is(get_se(), "NxtSE"), 
                    "Please build experiment first"))
            })

            req(is(get_se(), "NxtSE"))
            se <- isolate(get_se())
            
            # Slow step - populate genes
            req(get_ref())
            ref <- isolate(get_ref())
            withProgress(message = 'Populating gene list...', value = 0, {
                .server_cov_update_genes(session, ref$geneList)
            })
            
            settings_Cov$event.ranges <- as.data.table(
                coord2GR(rowData(se)$EventRegion))
            settings_Cov$event.ranges$EventName <- rowData(se)$EventName
            
            # refresh tracks and conditions here
            # - TODO

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
            
            .server_cov_update_events_list(session, res$EventName,
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
            chr_rd(), start_rd(), end_rd()
        ), {
            settings_Cov$localeTrigger <- runif(1)
        })

# Plot options triggers
        observeEvent(list(
            input$plot_Jn_cov, input$plot_key_iso, input$condense_cov,
        ), {
            settings_Cov$optionsTrigger <- runif(1)
        })

# Tracks selection trigger
        tracks_r <- reactive({
            req(input$track_table)
            hot_to_r(input$track_table) 
        })
        tracks_rd <- tracks_r %>% debounce(3000)
        diff_r <- reactive({
            req(is_valid(input$diffA))
            req(is_valid(input$diffB))
            return(c(input$diffA, input$diffB))
        })
        diff_rd <- diff_r %>% debounce(1000)

        observeEvent(list(
            tracks_rd(), diff_r(),
        ), {
            settings_Cov$tracksTrigger <- runif(1)
        })

# Start/End sanity check when these values are changed; also updates zoom
        observeEvent(chr_rd(), {
            seqInfo <- get_ref()$seqInfo[chr_rd()]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)
            output <- .server_cov2_check_start_end()
                input, session, output, seqmax)
        })
        observeEvent(list(start_rd(), end_rd()), {
            req(input$chr_cov, input$chr_cov %in% names(get_ref()$seqInfo))
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            seqmax <- as.numeric(GenomeInfoDb::seqlengths(seqInfo))
            req(seqmax > 50)
            output <- .server_cov2_check_start_end()
                input, session, output, seqmax)
        })

# Locate by Gene or Event triggers change in chr/start/end
        observeEvent(input$events_cov, {
            req(input$events_cov)
            req(input$events_cov != "(none)")
            
            events_id_view <- settings_Cov2$event.ranges[
                get("EventName") == input$events_cov]
            .server_cov_locate_events(input, session, events_id_view)
        })
        observeEvent(input$genes_cov, {
            req(input$genes_cov)
            req(input$genes_cov != "(none)")
            
            gene_id_view <- get_ref()$geneList[
                get("gene_display_name") == input$genes_cov]
            .server_cov_locate_genes(input, session, gene_id_view)
        })

# Sets start and end when zooming in / out
        observeEvent(input$zoom_out_cov, {
            req(input$zoom_out_cov, input$chr_cov) 
            req(input$chr_cov %in% names(get_ref()$seqInfo))
            
            seqInfo <- get_ref()$seqInfo[input$chr_cov]
            output <- .server_cov_zoom_out(
                input, output, session, seqInfo, settings_Cov)
        })
        observeEvent(input$zoom_in_cov, {
            req(input$zoom_in_cov)
            
            output <- .server_cov_zoom_in(input, output, session, settings_Cov)
        })

######################### Tracks management ####################################

# colData update changes condition drop-down box
        observeEvent(get_se(), {
            colData <- isolate(colData(get_se()))

            conditions <- colnames(colData)
            updateSelectInput(session = session, inputId = "condition_cov", 
                choices = c("(Individual Samples)", colnames(colData))
            )
        }

# change in condition will update tracks table
        observeEvent(input$condition_cov, {
            req(input$condition_cov)
            req(get_se())
            colData <- as.data.frame(isolate(colData(get_se())))

            if(input$condition_cov == "(Individual Samples)") {
                settings_Cov$trackTable <- data.frame(
                    names = colnames(get_se()),
                    id = 0
                )
            } else if(input$condition_cov %in% colnames(colData)) {
                condOptions <- unique(unname(unlist(
                    colData[, input$condition_cov])))
                settings_Cov$trackTable <- data.frame(
                    names = condOptions,
                    id = 0
                )
            }
        })

############################ Locale trigger ####################################

        observeEvent(settings_Cov$localeTrigger, {
            req(get_se())
            req(input$condition_cov)
            req(is_valid(settings_Cov$trackTable))

            args <- isolate(settings_Cov$dataObj)
            condName <- isolate(input$condition_cov)
            trackDF <- isolate(settings_Cov$trackTable)
            req("names" %in% colnames(trackDF))
            trackSamples <- trackDF$names
            
            if(!all(c("limit_start", "limit_end") %in% names(args))) {
                # Likely not a valid cDO, regenerate it
                settings_Cov$dataObj <- getCoverageData(
                    isolate(get_se()),
                    seqname = input$chr_cov,
                    start = input$start_cov,
                    end = input$end_cov,
                    strand = input$strand_cov,
                    tracks = colnames(get_se())
                )
            } else {
                # Check validity of chr/start/end call
                

            }
            
        })
        
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

