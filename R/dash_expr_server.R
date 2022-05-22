server_expr <- function(
        id, refresh_tab, volumes, get_threads_reactive, 
		get_memmode_reactive,
		limited = FALSE
) {
    moduleServer(id, function(input, output, session) {
        ns = NS(id)
        
        # Instantiate settings
        settings_expr <- setreactive_expr()
        
        # Updates to NxtSE loaded object. Refresh tab does nothing special
        observeEvent(list(refresh_tab(), settings_expr$se), {
            output <- .server_expr_parse_collate_path(limited = limited,
                settings_expr = reactiveValuesToList(settings_expr), 
                output = output)
        })
        
        # Directory and file handling
        observe({
            shinyDirChoose(input, "dir_reference_path_load", 
                roots = volumes(), session = session)
            shinyDirChoose(input, "dir_bam_path_load", 
                roots = volumes(), session = session)
            shinyDirChoose(input, "dir_sw_path_load", 
                roots = volumes(), session = session)
            shinyDirChoose(input, "dir_collate_path_load", 
                roots = volumes(), session = session)
            shinyFileChoose(input, "file_expr_anno_load", 
                roots = volumes(), session = session)
            shinyFileChoose(input, "loadNxtSE_RDS", 
                roots = volumes(), session = session)
        })
        observeEvent(input$dir_reference_path_load, {
            req(input$dir_reference_path_load)
            settings_expr$ref_path = parseDirPath(volumes(), 
                input$dir_reference_path_load)
        })
        observeEvent(input$dir_bam_path_load, {
            req(input$dir_bam_path_load)
            settings_expr$bam_path = parseDirPath(volumes(), 
                input$dir_bam_path_load)
        })
        observeEvent(input$dir_sw_path_load, {
            req(input$dir_sw_path_load)
            settings_expr$sw_path = parseDirPath(volumes(), 
                input$dir_sw_path_load)
        })
        observeEvent(input$file_expr_anno_load, {
            req(input$file_expr_anno_load)
            settings_expr$anno_file = as.character(
                parseFilePaths(volumes(), input$file_expr_anno_load))
        })
        observeEvent(input$loadNxtSE_RDS, {
            req(input$loadNxtSE_RDS)
            RDSfile <- as.character(
                parseFilePaths(volumes(), input$loadNxtSE_RDS))
            tmpse <- readRDS(RDSfile)
            if(!is(tmpse, "NxtSE")) {
                .load_NxtSE_sweetalert_finish(session)
                settings_expr$se <- tmpse
            } else {
                .load_NxtSE_sweetalert_error(session)
            }
            rm(tmpse)
        })
        observeEvent(input$dir_collate_path_load, {
            req(input$dir_collate_path_load)
            settings_expr$collate_path = parseDirPath(volumes(), 
                input$dir_collate_path_load)
        })
        
        # Experiment I/O - sync between files and anno
        observeEvent(settings_expr$df.files, {
            req(settings_expr$df.files)
            req(is(settings_expr$df.files, "data.frame"))
            req("sample" %in% colnames(settings_expr$df.files))
            settings_expr$df.anno <- .server_expr_sync_df(
                settings_expr$df.files, settings_expr$df.anno)
        })
        observeEvent(settings_expr$df.anno, {
            req(settings_expr$df.anno)
            req(is(settings_expr$df.anno, "data.frame"))
            req("sample" %in% colnames(settings_expr$df.anno))
            req(settings_expr$df.files)
            settings_expr$df.files <- .server_expr_sync_df(
                settings_expr$df.anno, settings_expr$df.files)
            # If annotations are added, validate NxtSE Object
            output <- .server_expr_parse_collate_path(
                limited = limited,
                settings_expr = reactiveValuesToList(settings_expr), 
                output = output
            )
        })
        
        # Experiment I/O - sync between user input and data frames
        observeEvent(input$hot_files_expr,{
            req(input$hot_files_expr)
            settings_expr$df.files = hot_to_r(input$hot_files_expr) 
        })
        observeEvent(input$hot_anno_expr,{
            req(input$hot_anno_expr)
            settings_expr$df.anno = hot_to_r(input$hot_anno_expr)
        })
        output$hot_files_expr <- renderRHandsontable({
            .server_expr_gen_HOT(settings_expr$df.files, enable_select = TRUE)
        })
        output$hot_anno_expr <- renderRHandsontable({
            .server_expr_gen_HOT(settings_expr$df.anno)
        })
        
        # Experiment I/O - saves and loads to NxtSE project directory
        observeEvent(input$save_expr,{
            req(input$save_expr)
            .server_expr_save_expr(reactiveValuesToList(settings_expr), session)
            settings_expr$df.files_savestate <- settings_expr$df.files
            settings_expr$df.anno_savestate <- settings_expr$df.anno
            # Validate NxtSE Object
            output <- .server_expr_parse_collate_path(
                limited = limited,
                settings_expr = reactiveValuesToList(settings_expr), 
                output = output
            )
        })
        
        # Click "Load Experiment" after setting collateData output path
        observeEvent(input$load_expr,{
            req(input$load_expr)
            output <- .server_expr_load_expr(
                reactiveValuesToList(settings_expr), session, output)
            # Update status boxes
            output <- .server_expr_parse_collate_path(
                limited = limited,
                settings_expr = reactiveValuesToList(settings_expr), 
                output = output
            )
        })

        # Toggle to Annotations
        observeEvent(input$add_anno, {
            updateRadioGroupButtons(session, inputId = "hot_switch_expr", 
                selected = "Annotations")
        })
        output$newcol_expr <- renderUI({
            textInput(ns("newcolumnname_expr"), "New Column Name (to add)", 
                sprintf("newcol%s", 1 + ncol(settings_expr$df.anno))
            )
        })
        
        # Add new annotation column
        observeEvent(input$addcolumn_expr, {
            req(input$addcolumn_expr)
            df <- isolate(settings_expr$df.anno)
            newcolumn <- eval(parse(text=sprintf('%s(nrow(df))', 
                isolate(input$type_newcol_expr))))
            settings_expr$df.anno <- data.table::setnames(
                cbind(df, newcolumn, stringsAsFactors=FALSE), 
                c(names(df), isolate(input$newcolumnname_expr))
            )
        })
        
        # Remove annotation column
        observeEvent(input$removecolumn_expr, {
            req(input$removecolumn_expr)
            DT <- as.data.table(isolate(settings_expr$df.anno))
            if(isolate(input$newcolumnname_expr) %in% colnames(DT)) {
                message("removing column")
                DT[, c(input$newcolumnname_expr) := NULL]
                settings_expr$df.anno = as.data.frame(DT)
            }
        })
        
        # Clearing Selections
        observeEvent(input$dir_collate_path_clear, {
            req(input$dir_collate_path_clear)
            settings_expr$collate_path = ""
        })
        observeEvent(input$clearLoadRef,{
            settings_expr$ref_path = ""
            output <- .server_expr_clear_ref(output)   
        })
        observeEvent(input$clear_expr, {
			settings_expr$ref_path = ""
            settings_expr$bam_path = ""
            settings_expr$sw_path = ""
            settings_expr$anno_file = ""
            settings_expr$collate_path = ""
            settings_expr$df.files = c()
            settings_expr$df.anno = c()
            settings_expr$se = NULL
        })    
        
        # Event when reference directory is set
        observeEvent(settings_expr$ref_path, {
            path_selected           <- is_valid(settings_expr$ref_path)
            path                    <- settings_expr$ref_path
            settings_expr$ref_path  <- .server_expr_check_ref_path(
                                        settings_expr$ref_path)
            output                  <- .server_expr_parse_ref_path(
                                        settings_expr$ref_path, output)
            
            if(is_valid(settings_expr$ref_path)) {
                settings_expr$ref_settings <- readRDS(
                    file.path(settings_expr$ref_path, "settings.Rds"))
                .server_expr_ref_load_success(session, settings_expr$ref_path)
            } else if(path_selected) {
                settings_expr$ref_settings <- NULL
                .server_expr_ref_load_fail(session, path)      
            }
        })

        # Event when BAM directory is set
        observeEvent(settings_expr$bam_path,{
            settings_expr$df.files <- Expr_Load_BAMs(
                settings_expr$df.files, settings_expr$bam_path, session)
            output$bam_expr_infobox <- Expr_BAM_update_status(
                settings_expr$df.files, settings_expr$bam_path,
                settings_expr$collate_path)
            output$txt_bam_path_load <- renderText(
                settings_expr$bam_path)
        })

        # Event when processBAM output directory is set
        observeEvent(settings_expr$sw_path,{
            settings_expr$df.files <- Expr_Load_SW(
                settings_expr$df.files, settings_expr$sw_path)
            output <- .server_expr_check_sw_path(settings_expr$df.files, 
                settings_expr$sw_path, output)
            output$txt_sw_path_expr <- renderText(
                settings_expr$sw_path)
        })

        # Event when Annotation file is set
        observeEvent(settings_expr$anno_file,{
            req(settings_expr$anno_file)
            settings_expr$df.anno <- Expr_Load_Anno(settings_expr$df.anno,
                settings_expr$df.files, settings_expr$anno_files)
        })

        # Event when NxtSE output directory is set
        observeEvent(settings_expr$collate_path, {
            if(
                    is_valid(settings_expr$collate_path) && 
                    file.exists(
                        file.path(settings_expr$collate_path, "colData.Rds")
                    )
            ) {
                colData.Rds = readRDS(file.path(
                    settings_expr$collate_path, "colData.Rds"))
                if(all(c("df.anno", "df.files") %in% names(colData.Rds))) {
                    settings_expr$df.files              <- colData.Rds$df.files
                    settings_expr$df.anno               <- colData.Rds$df.anno
                    settings_expr$df.files_savestate <- settings_expr$df.files
                    settings_expr$df.anno_savestate <- settings_expr$df.anno
                    if("bam_path" %in% names(colData.Rds)) {
                        settings_expr$bam_path <- colData.Rds$bam_path
                    }
                    if("sw_path" %in% names(colData.Rds)) {
                        settings_expr$sw_path <- colData.Rds$sw_path
                    }
                }
            }
            output <- .server_expr_parse_collate_path(
                limited = limited,
                settings_expr = reactiveValuesToList(settings_expr), 
                output = output
            )
            output$txt_NxtSE_path_load <- renderText(
                settings_expr$collate_path)
        })

        # Running processBAM
        observeEvent(input$run_pb_expr,{
            req(input$run_pb_expr)
            settings_expr$selected_rows <- Expr_PB_initiate_run(
                input, session, 
                get_threads_reactive(), 
                isolate(reactiveValuesToList(settings_expr))
            )
        })
        observeEvent(input$pb_confirm, {
            if(input$pb_confirm == FALSE) {
                settings_expr$selected_rows = c()
                return()
            } else {
                Expr_PB_actually_run(
                    input, session, get_threads_reactive(), 
                    isolate(reactiveValuesToList(settings_expr))
                )
            }
            settings_expr$selected_rows = c()
            settings_expr$df.files <- Expr_Load_SW(
                settings_expr$df.files, settings_expr$sw_path)
            output <- .server_expr_check_sw_path(settings_expr$df.files, 
                settings_expr$sw_path, output)
            output <- .server_expr_parse_collate_path(limited = limited,
                settings_expr = reactiveValuesToList(settings_expr), 
                output = output)
        })

        # Running collateData
        observeEvent(input$run_collate_expr, {
            req(input$run_collate_expr)
            req(settings_expr$df.files)
            Experiment = na.omit(as.data.table(
                settings_expr$df.files[, c("sample", "sw_file", "cov_file")]
            ))
            reference_path = settings_expr$ref_path
            output_path = settings_expr$collate_path
            if(Expr_collateData_Validate_Vars(
                    session, Experiment, reference_path, output_path
            )) {
                withProgress(
                        message = 'Collating SpliceWiz (processBAM) output', 
                        value = 0, 
                {
                    collateData(
                        Experiment, reference_path, output_path, 
                        n_threads = get_threads_reactive(),
						lowMemoryMode = get_memmode_reactive()
                    )
                })
                Expr_Update_colData(
                    settings_expr$collate_path, 
                    settings_expr$df.anno, settings_expr$df.files, 
                    settings_expr$bam_path, settings_expr$sw_path, 
                    session, post_collateData = TRUE
                )   # saves / updates expr
                output <- .server_expr_parse_collate_path(
                    limited = limited,
                    settings_expr = reactiveValuesToList(settings_expr), 
                    output = output
                )    # updates status boxes
            }
        })

        # Running makeSE (Only available on limited == TRUE)
        observeEvent(input$build_expr, {
            if(
                    is_valid(settings_expr$collate_path) &&
                    file.exists(file.path(
                        settings_expr$collate_path, "colData.Rds"))
            ) {
                colData = as.data.table(settings_expr$df.anno)
                withProgress(message = 'Loading NxtSE object', value = 0, {
                    tryCatch({
                        settings_expr$se = makeSE(
                            settings_expr$collate_path, colData,
                            realize = TRUE
                        )
                        .makeSE_sweetalert_finish(session)
                    }, error = function(e) {
                        .makeSE_sweetalert_error(session)
                    })
                })
            }
        })

        observeEvent(input$saveNxtSE_RDS, {    
            req(settings_expr$se)
            selectedfile <- parseSavePath(volumes(), input$saveNxtSE_RDS)
            req(selectedfile$datapath)
            saveRDS(settings_expr$se, selectedfile$datapath)
        })

    # End of Server function
        return(settings_expr)
    })
}

# Clear all info boxes
.server_expr_clear_ref <- function(output) {
    # output$fasta_source_infobox <- renderInfoBox(infoBox(""))
    # output$gtf_source_infobox <- renderInfoBox(infoBox(""))
    # output$mappa_source_infobox <- renderInfoBox(infoBox(""))
    # output$NPA_source_infobox <- renderInfoBox(infoBox(""))
    # output$BL_source_infobox <- renderInfoBox(infoBox(""))
    output$txt_reference_path_load <- renderText("")
    output$ref_expr_infobox <- renderUI(ui_infobox_ref(""))
    return(output)
}

# Check path contains valid SpliceWiz reference
.server_expr_check_ref_path <- function(ref_path) {
    if(is_valid(ref_path)) {
        ref_settings_file <- file.path(ref_path, "settings.Rds")
        if(file.exists(ref_settings_file)) {
            ref_settings = readRDS(ref_settings_file)
            if("reference_path" %in% names(ref_settings)) {
                return(ref_path)
            }
        }
    }
    return("")
}

# Register ref_path into server
.server_expr_parse_ref_path <- function(ref_path, output) {
    if(is_valid(ref_path)) {
        ref_settings_file <- file.path(ref_path, "settings.Rds")
        ref_settings <- readRDS(ref_settings_file)
        output <- .server_expr_load_ref(ref_settings, 
            output)
        output$txt_reference_path_load <- renderText(
            ref_path)
        output$ref_expr_infobox <- renderUI(ui_infobox_ref(
            ref_settings_file))
    } else {
        output <- .server_expr_clear_ref(output)
    }
    return(output)
}

# Filter df2 by the samples in df1
.server_expr_sync_df <- function(df1, df2) {
    if(!is_valid(df2)) {
        return(data.frame(sample = df1$sample, stringsAsFactors = FALSE))
    } else {
        DT1 = as.data.table(df1)
        DT2 = as.data.table(df2)
        samples = DT1[, "sample"]
        return(as.data.frame(DT2[samples, on = "sample"]))
    }
}

# Generate rHOT from df
.server_expr_gen_HOT <- function(df, enable_select = FALSE) {
    if(is_valid(df) && is(df, "data.frame")) {
        rhandsontable(df, useTypes = TRUE, stretchH = "all",
            selectCallback = enable_select)
    } else {
        NULL
    }
}

# Load settings.Rds from SpliceWiz reference to populate status boxes
.server_expr_load_ref <- function(ref_settings, output) {
    ah <- ah_genome_record <- ah_gtf_record <- NULL
    fasta <- gtf <- mappa <- nonPA <- Black <- NULL
    if(
            "ah_genome" %in% names(ref_settings) &&
            is_valid(ref_settings[["ah_genome"]])
    ) {
        ah <- AnnotationHub()
        ah_genome_record <- tryCatch({
            basename(ah$sourceurl[
                which(names(ah) == ref_settings[["ah_genome"]])])
        }, error = function(e) NULL)
    }
    if(
            "ah_transcriptome" %in% names(ref_settings) &&
            is_valid(ref_settings[["ah_transcriptome"]])
    ) {
        if(is.null(ah)) ah <- AnnotationHub()
        ah_gtf_record <- tryCatch({
            basename(ah$sourceurl[
                which(names(ah) == ref_settings[["ah_transcriptome"]])])
        }, error = function(e) NULL)
    }
    if(is.null(ah_genome_record) && "fasta_file" %in% names(ref_settings)) {
        fasta <- basename(ref_settings[["fasta_file"]])
    }
    if(is.null(ah_gtf_record) && "gtf_file" %in% names(ref_settings)) {
        gtf <- basename(ref_settings[["gtf_file"]])
    }
    if("MappabilityRef" %in% names(ref_settings)) {
        mappa <- basename(ref_settings[["MappabilityRef"]])
    }
    if("nonPolyARef" %in% names(ref_settings)) {
        nonPA <- basename(ref_settings[["nonPolyARef"]])
    }
    if("BlacklistRef" %in% names(ref_settings)) {
        Black <- basename(ref_settings[["BlacklistRef"]])
    }   
    output <- .server_expr_load_ref_genome(output, ah_genome_record, fasta)
    output <- .server_expr_load_ref_gtf(output, ah_gtf_record, gtf)
    output <- .server_expr_load_ref_misc(output, mappa, nonPA, Black)
    return(output)
}

# Add genome FASTA profile into server
.server_expr_load_ref_genome <- function(output, ah_genome_record, fasta) {
    if(is_valid(ah_genome_record)) {
        output$fasta_source_infobox <- renderInfoBox({
            infoBox("Genome - AnnotationHub", "", ah_genome_record,
                icon = icon("dna", lib = "font-awesome"), 
                color = "green")
        })      
    } else if(is_valid(fasta)) {
        output$fasta_source_infobox <- renderInfoBox({
            infoBox("Genome - User FASTA", "",
                fasta, 
                icon = icon("dna", lib = "font-awesome"),
                color = "green")
        })
    } else {
        output$fasta_source_infobox <- renderInfoBox({
            infoBox("Genome - INVALID", "",
                "", 
                icon = icon("dna", lib = "font-awesome"),
                color = "red")
        })
    }
    return(output)
}

# Add annotation GTF profile into server
.server_expr_load_ref_gtf <- function(output, ah_gtf_record, gtf) {
    if(is_valid(ah_gtf_record)) {
        output$gtf_source_infobox <- renderInfoBox({
            infoBox("Gene Annotation - AnnotationHub",  "", ah_gtf_record,
                icon = icon("book-medical", lib = "font-awesome"),
                color = "orange")
        })               
    } else if(is_valid(gtf)) {
        output$gtf_source_infobox <- renderInfoBox({
            infoBox("Gene Annotation - User GTF",  "",
                gtf, 
                icon = icon("book-medical", lib = "font-awesome"),
                color = "orange")
        })
    } else {
        output$gtf_source_infobox <- renderInfoBox({
            infoBox("Gene Annotation - INVALID", "",
                "", 
                icon = icon("book-medical", lib = "font-awesome"),
                color = "red")
        })
    }
    return(output)
}

# Add miscellaneous reference data into server
.server_expr_load_ref_misc <- function(output, mappa, nonPA, Black) {
    if(is_valid(mappa)) {
        output$mappa_source_infobox <- renderInfoBox({
            infoBox("Mappability", "",
                mappa, 
                icon = icon("map", lib = "font-awesome"),
                color = "blue")
        })
    } else {
        output$mappa_source_infobox <- renderInfoBox({
            infoBox("Mappability", "",
                "NOT USED", 
                icon = icon("map", lib = "font-awesome"),
                color = "blue")
        })  
    }
    if(is_valid(nonPA)) {
        output$NPA_source_infobox <- renderInfoBox({
            infoBox("Non-PolyA", "",
                nonPA, 
                icon = icon("font", lib = "font-awesome"),
                color = "purple")
        })
    } else {
        output$NPA_source_infobox <- renderInfoBox({
            infoBox("Non-PolyA", "",
                "NOT USED", 
                icon = icon("font", lib = "font-awesome"),
                color = "purple")
        })
    }
    if(is_valid(Black)) {
        output$BL_source_infobox <- renderInfoBox({
            infoBox("BlackList", "",
                Black, 
                icon = icon("list-alt", lib = "font-awesome"),
                color = "red")
        })
    } else {
        output$BL_source_infobox <- renderInfoBox({
            infoBox("BlackList", "",
                "NOT USED", 
                icon = icon("list-alt", lib = "font-awesome"),
                color = "red")
        })
    }
    return(output)
}

# Given a BAM path, load BAM files to populate experiment
Expr_Load_BAMs <- function(df.files, bam_path, session) {
    if(!is_valid(bam_path)) return(df.files)

    # First assume bams are named by subdirectory names
    temp.DT <- findSamples(bam_path, suffix = ".bam", level = 1)
    if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
        temp.DT <- as.data.table(temp.DT)
        if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Assume subdirectory names designate sample names
        } else {
            temp.DT <- as.data.table(findSamples(
                bam_path, suffix = ".bam", level = 0))
            if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Else assume bam names designate sample names
            } else {
                sendSweetAlert(session = session,
                    title = "Incompatible BAM file names",
                    text = paste("Could not determine sample names.",
                        "Please ensure either BAMs are uniquely named by",
                        "sample name,",
                        "or its parent directories are uniquely named."
                    ), type = "error")
                temp.DT <- NULL
            }
        }
    } else {
        sendSweetAlert(session = session, 
            title = "No BAM files found",
            text = "No BAM files found", type = "error")            
        temp.DT = NULL
    }
    # compile experiment df with bam paths
    if(!is.null(temp.DT) && nrow(temp.DT) > 0)  {
        colnames(temp.DT)[2] <- "bam_file"
        if(is_valid(df.files)) {
            df.files <- update_data_frame(df.files, temp.DT)
        } else {
            DT <- data.table(sample = temp.DT$sample,
                bam_file = "", sw_file = "", cov_file = "")
            DT[temp.DT, on = "sample", c("bam_file") := get("i.bam_file")]
            df.files <- as.data.frame(DT)
        }
        return(df.files)
    } else {
        return(df.files)
    }
}

Expr_BAM_update_status <- function(df.files, bam_path, collate_path) {
    if(is_valid(df.files)) {
        if(
                is_valid(bam_path) &&
                "bam_file" %in% colnames(df.files) && 
                all(file.exists(df.files$bam_file))
        ) {
            return(renderUI(ui_infobox_bam(bam_path, df.files$bam_file)))
        } else if(
                "sw_file" %in% colnames(df.files) && 
                all(file.exists(df.files$sw_file))
        ) {
            return(renderUI(ui_infobox_bam(bam_path, escape = TRUE)))
        } else if(
                is_valid(collate_path) && 
                file.exists(file.path(collate_path, "colData.Rds"))
        ) {
            return(renderUI(ui_infobox_bam(bam_path, escape = TRUE)))
        } else if("bam_file" %in% colnames(df.files)) {
            return(renderUI(ui_infobox_bam(bam_path, df.files$bam_file)))
        }
    } else {
        return(renderUI(ui_infobox_bam(bam_path)))        
    } 
}

Expr_Load_SW <- function(df.files, sw_path) {
    if(!is_valid(sw_path)) return(df.files)
    # merge splicewiz paths
    temp.DT <- findSamples(sw_path, suffix = ".txt.gz", level = 0)
    if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
        temp.DT <- as.data.table(temp.DT)
        if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Assume output names designate sample names
        } else {
            temp.DT <- as.data.table(findSamples(
                sw_path, suffix = ".txt.gz", level = 1))
            if(length(unique(temp.DT$sample)) == nrow(temp.DT)) {
            # Else assume subdirectory names designate sample names
            } else {
                temp.DT <- NULL
            }
        }
    } else {
        temp.DT <- NULL
    }
    if(!is.null(temp.DT) && nrow(temp.DT) > 0) {
        colnames(temp.DT)[2] <- "sw_file"
        if(is_valid(df.files)) {
            df.files <- update_data_frame(df.files, temp.DT)
        } else {
            DT <- data.table(sample = temp.DT$sample,
                bam_file = "", sw_file = "", cov_file = "")
            DT[temp.DT, on = "sample", c("sw_file") := get("i.sw_file")] 
            df.files <- as.data.frame(DT)      
        }   
    }
    temp.DT2 <- findSamples(sw_path, suffix = ".cov", level = 0)
    if(!is.null(temp.DT2) && nrow(temp.DT2) > 0) {
        temp.DT2 <- as.data.table(temp.DT2)
        if(length(unique(temp.DT2$sample)) == nrow(temp.DT2)) {
            # Assume output names designate sample names
        } else {
            temp.DT2 <- as.data.table(findSamples(
                sw_path, suffix = ".cov", level = 1))
            if(length(unique(temp.DT2$sample)) == nrow(temp.DT2)) {
        # Else assume subdirectory names designate sample names
            } else {
                temp.DT2 <- NULL
            }
        }
    } else {
        temp.DT2 <- NULL
    }
# compile experiment df with splicewiz paths
    if(!is.null(temp.DT2) && nrow(temp.DT2) > 0) {
        colnames(temp.DT2)[2] <- "cov_file"
        df.files <- update_data_frame(df.files, temp.DT2)
    }
    return(df.files)
}

# Brings a prompt message asking do you really want to run processBAM
Expr_PB_initiate_run <- function(input, session, n_threads, settings_expr) {
    if(!is_valid(settings_expr$df.files)) {
        sendSweetAlert(session = session, type = "error",
            title = "No bam files in experiment",
            text = "Please select bam folder and select bam files")
        return()
    }
    if(!("bam_file" %in% colnames(settings_expr$df.files))) {
        sendSweetAlert(session = session, type = "error",
            title = "No bam files in experiment",
            text = "Please select bam folder and select bam files")
        return()
    }
    if(!is_valid(input$hot_files_expr_select$select$r)) {
        sendSweetAlert(session = session, type = "error",
            title = "No BAM files selected", 
            text = "Please highlight cells of bam files to run processBAM")
        return()        
    }
    selected_rows <- seq(input$hot_files_expr_select$select$r,
        input$hot_files_expr_select$select$r2)
    selected_cols <- seq(input$hot_files_expr_select$select$c,
        input$hot_files_expr_select$select$c2)
    bam_col <- which(colnames(settings_expr$df.files) == "bam_file")
    bam_files <- settings_expr$df.files$bam_file[selected_rows]
    if(!is_valid(settings_expr$ref_path)) {
        sendSweetAlert(session = session,
            title = "Missing Reference", type = "error",
            text = "Please load Reference before running processBAM")
    } else if(!(bam_col %in% selected_cols)) {
        sendSweetAlert(session = session, type = "error",
            title = "No BAM files selected",
            text = "Please highlight cells of bam files to run processBAM")
    } else if(!all(file.exists(bam_files))) {
        sendSweetAlert(session = session,
            title = "Missing BAMs", type = "error",
            text = "Please check all selected bam files exist")          
    } else if(!file.exists(file.path(
            settings_expr$ref_path, "SpliceWiz.ref.gz"))) {
        sendSweetAlert(session = session, type = "error",
            title = "Missing SpliceWiz Reference",
            text = "SpliceWiz.ref.gz is missing")
    } else if(
            !is_valid(settings_expr$sw_path) || 
            !dir.exists(settings_expr$sw_path)
    ) {
        sendSweetAlert(session = session, type = "error",
            title = "Missing SpliceWiz (processBAM) output path",
            text = "Please set SpliceWiz (processBAM) output path")
    } else {
        # n_threads <- min(n_threads, length(selected_rows))
        # if(n_threads < length(selected_rows)) {
            # n_rounds <- ceiling(length(selected_rows) / n_threads)
            # n_threads <- ceiling(length(selected_rows) / n_rounds)
        # }
        msg <- paste("Run processBAM on", length(selected_rows), "samples",
            # "Estimated runtime", 
                # 10 * ceiling(length(selected_rows) / n_threads), "minutes"
            "using", n_threads, "threads"#,
            # "threads (10min per BAM @ 100 million reads per sample)"
        )
        ask_confirmation(inputId = "pb_confirm", type = "warning", 
            title = msg, btn_labels = c("Cancel", "Run processBAM"),
            btn_colors = c("#00BFFF", "#FE2E2E"))
        return(selected_rows)
    }
}

# After user confirms, actually call processBAM
Expr_PB_actually_run <- function(input, session, n_threads, settings_expr) {

    withProgress(message = 'Running SpliceWiz (processBAM)', value = 0, {
        i_done <- 0
        incProgress(0.001, 
            message = paste('Running SpliceWiz (processBAM)',
                i_done, "of", length(settings_expr$selected_rows), "done")
        )
        for(i in settings_expr$selected_rows) {
            processBAM(
                bamfiles = settings_expr$df.files$bam_file[i],
                sample_names = settings_expr$df.files$sample[i],
                reference_path = settings_expr$ref_path,
                output_path = settings_expr$sw_path,
                n_threads = n_threads,
                run_featureCounts = FALSE,
                verbose = TRUE                    
            )
            i_done <- i_done + 1
            incProgress(1 / length(settings_expr$selected_rows), 
                message = paste(i_done, "of", 
                    length(settings_expr$selected_rows), "done")
            )
        }
    })

    sendSweetAlert(
        session = session,
        title = "SpliceWiz (processBAM) run completed",
        type = "success"
    )
}

# Check SpliceWiz path contains SpliceWiz output or not
.server_expr_check_sw_path <- function(df.files, sw_path, output) {
    if(is_valid(df.files) && "sw_file" %in% colnames(df.files)) {
        sw_files <- df.files$sw_file
    } else {
        sw_files <- NULL
    }
    output$pb_expr_infobox <- renderUI({
        ui_infobox_pb(sw_path, sw_files)
    })
    return(output)
}

# Load annotation file
Expr_Load_Anno <- function(df.anno, df.files, anno_file, session) {
    temp.DT <- tryCatch(fread(anno_file), error = function(e) NULL)
    if(!is_valid(temp.DT)) return(df.anno)
    if(nrow(temp.DT) == 0) return(df.anno)
    if(!("sample" %in% colnames(temp.DT))) {
        sendSweetAlert(
            session = session,
            title = "Error in Annotation file",
            text = "'sample' must be the name of the first column",
            type = "error"
        )
        return(df.anno)
    }
    
    files_header <- c("bam_file", "sw_file", "cov_file")
    anno_header <- names(temp.DT)[!(names(temp.DT) %in% files_header)]
    temp.DT.files <- copy(temp.DT)
    if(length(anno_header) > 0) temp.DT.files[, c(anno_header) := NULL]
    if(is_valid(df.files)) {
        df.files <- update_data_frame(df.files, temp.DT.files)
    } else {
        DT <- data.table(
            sample = temp.DT$sample, bam_file = "", 
            sw_file = "", cov_file = ""
        )
        df.files <- update_data_frame(DT, temp.DT.files)
    }
    temp.DT.anno <- copy(temp.DT)
    files_header_exist <- intersect(files_header, names(temp.DT))
    if(length(files_header_exist) > 0) {
        temp.DT.anno[, c(files_header_exist):= NULL]
    }
    if(is_valid(df.anno)) {
        df.anno <- update_data_frame(df.anno, temp.DT.anno)
    } else {
        df.anno <- temp.DT.files
    }
    return(df.anno)
}

# Check if savestate df is identical to loaded df
.server_expr_check_savestate <- function(settings_expr) {
    return(
        identical(settings_expr$df.anno_savestate, settings_expr$df.anno) &&
        identical(settings_expr$df.files_savestate, settings_expr$df.files)
    )
}

.server_expr_parse_collate_path <- function(limited, settings_expr, output) {
    if(limited) {
        return(.server_expr_parse_collate_path_limited(settings_expr, output))
    } else {
        return(.server_expr_parse_collate_path_full(settings_expr, output))
    }
}

# Checks collate path and report status
.server_expr_parse_collate_path_limited <- function(settings_expr, output) {
    if(is_valid(settings_expr$se)) {
        if(
                ncol(settings_expr$df.anno) > 1 &&
                .server_expr_check_savestate(settings_expr)
        ) {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(3, "NxtSE Loaded"))
        } else if(ncol(settings_expr$df.anno) > 1) {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(2, "NxtSE Loaded",
                    "Don't forget to save your experiment"))
        } else {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(1, "NxtSE Loaded",
                    "Consider adding one or more conditions to Annotations"))
        }
    } else if(
            is_valid(settings_expr$collate_path) &&
            file.exists(file.path(
                settings_expr$collate_path, "NxtSE.rds"))
    ) {
        if(
                ncol(settings_expr$df.anno) > 1 && 
                .server_expr_check_savestate(settings_expr)
        ) {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(2, "NxtSE ready to load",
                    "Click `Build SummarizedExperiment`"))
        } else if(ncol(settings_expr$df.anno) > 1) {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(1, "NxtSE ready to load",
                    # "Click `Build SummarizedExperiment`",
                    "Don't forget to save your experiment"))
        } else {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(1, "NxtSE ready to load",
                    "Consider adding conditions to Annotations"))
        }
    } else if(
            is_valid(settings_expr$collate_path) &&
            is_valid(settings_expr$df.files) &&
            all(file.exists(settings_expr$df.files$sw_file))
    ) {
        output$se_expr_infobox <- renderUI(
            ui_infobox_expr(1, "NxtSE not collated",
                "Run collateData via Experiment tab"))
    } else if(is_valid(settings_expr$collate_path)) {
        output$se_expr_infobox <- renderUI(
            ui_infobox_expr(0,
            submsg = "Run processBAM and collateData via the Experiment tab"))
    } else {
        output$se_expr_infobox <- renderUI(
            ui_infobox_expr(0,
            submsg = "Select output directory of collated data"))
    }
    return(output)
}

# Checks collate path and report status
.server_expr_parse_collate_path_full <- function(settings_expr, output) {
    if(
            is_valid(settings_expr$collate_path) &&
            file.exists(file.path(settings_expr$collate_path, "NxtSE.rds"))
    ) {
        if(.server_expr_check_savestate(settings_expr)) {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(3, "NxtSE ready to load", 
                    "Load via Analysis -> Load Experiment"))
        } else {
            output$se_expr_infobox <- renderUI(
                ui_infobox_expr(2, "NxtSE ready to load", 
                    "Don't forget to save your experiment"))
        }
    } else if(
            is_valid(settings_expr$collate_path) &&
            is_valid(settings_expr$df.files) &&
            all(file.exists(settings_expr$df.files$sw_file))
    ) {
        output$se_expr_infobox <- renderUI(
            ui_infobox_expr(2, "Ready to collate experiment"))
    } else if(
            is_valid(settings_expr$collate_path) && 
            is_valid(settings_expr$df.files)
    ) {
        output$se_expr_infobox <- renderUI(
            ui_infobox_expr(1, "processBAM output files incomplete"))
    } else if(is_valid(settings_expr$collate_path)) {
        output$se_expr_infobox <- renderUI(ui_infobox_expr(0,
            paste("Selected path:", settings_expr$collate_path)))
    } else {
        output$se_expr_infobox <- renderUI(ui_infobox_expr(0,
            "Select path for NxtSE output"))
    }
    return(output)
}

# Save annotations to colData.Rds
.server_expr_save_expr <- function(settings_expr, session) {
    if(
        is_valid(settings_expr$collate_path)
    ) {
        colData.Rds = list(
            df.anno = settings_expr$df.anno,
            df.files = settings_expr$df.files,
            bam_path = settings_expr$bam_path,
            sw_path = settings_expr$sw_path
        )
        saveRDS(colData.Rds, file.path(settings_expr$collate_path, 
            "colData.Rds"))
        sendSweetAlert(
            session = session,
            title = paste("Annotations saved to", settings_expr$collate_path),
            type = "success"
        )
    } else {
        sendSweetAlert(
            session = session,
            title = "Annotations not saved; run collateData first!",
            type = "error"
        )
    }
}

.server_expr_load_expr <- function(settings_expr, session, output) {
    if(
            is_valid(settings_expr$collate_path) &&
            file.exists(
                file.path(settings_expr$collate_path, "colData.Rds")
            )
    ) {
        colData.Rds = readRDS(
            file.path(settings_expr$collate_path, "colData.Rds"))
        req_columns = c("df.anno", "df.files")
        if(all(req_columns %in% names(colData.Rds))) {
            settings_expr$df.files <- colData.Rds$df.files
            settings_expr$df.files_savestate <- settings_expr$df.files
            settings_expr$df.anno <- colData.Rds$df.anno
            settings_expr$df.anno_savestate <- settings_expr$df.anno
            if("bam_path" %in% names(colData.Rds)) {
                settings_expr$bam_path = colData.Rds$bam_path
            }
            if("sw_path" %in% names(colData.Rds)) {
                settings_expr$sw_path = colData.Rds$sw_path
            }
            output <- .server_expr_parse_collate_path(
                limited = limited,
                settings_expr = reactiveValuesToList(settings_expr), 
                output = output
            )
            .server_expr_load_alert_success(session, 
                settings_expr$collate_path)
        } else {
            .server_expr_load_alert_fail(session, 
                settings_expr$collate_path)
        }
    } else {
        .server_expr_load_alert_fail(session, 
            settings_expr$collate_path)
    }
    return(output)
}

# Check paths are legit before running collateData()
Expr_collateData_Validate_Vars <- function(
        session, Experiment, reference_path, output_path
) {
    if(!is_valid(reference_path)) {
        sendSweetAlert(
            session = session,
            title = "Missing Reference",
            text = "Please load Reference before running collateData",
            type = "error"
        )
        return(FALSE)
    } else if(!is_valid(output_path)) {
        sendSweetAlert(
            session = session,
            title = "Missing NxtSE Path",
            text = paste("Please select NxtSE path before",
                "running collateData"),
            type = "error"
        )
        return(FALSE)
    } else if(!dir.exists(output_path)) {
        sendSweetAlert(
            session = session,
            title = "Invalid NxtSE Path",
            text = "Please make sure NxtSE output path exists",
            type = "error"
        )
        return(FALSE)
    } else if(nrow(Experiment) == 0) {
        sendSweetAlert(
            session = session,
            title = "No samples found to collate Experiment",
            text = "Please load processBAM output of some samples",
            type = "error"
        )
        return(FALSE)
    }
    return(TRUE)
}

# Sends sweetAlerts to show whether collateData() has run successfully
Expr_Update_colData <- function(
        collate_path, df.anno, df.files, 
        bam_path, sw_path, session, 
        post_collateData = FALSE)
{
    if(file.exists(file.path(collate_path, "colData.Rds"))) {
        colData.Rds <- readRDS(file.path(collate_path, "colData.Rds"))
        if(all(colData.Rds$df.anno$sample %in% df.anno$sample)) {
            colData.Rds$df.anno <- df.anno
            colData.Rds$df.files <- df.files
            colData.Rds$bam_path <- bam_path
            colData.Rds$sw_path <- sw_path
            saveRDS(colData.Rds, file.path(collate_path, "colData.Rds"))
            if(post_collateData) {
                sendSweetAlert(
                    session = session,
                    title = "collateData run completed",
                    type = "success"
                )
            }
        } else {
            if(post_collateData) {
                sendSweetAlert(
                    session = session,
                    title = "collateData did not collate all samples",
                    type = "warning"
                )
            }
        }
    } else if(is_valid(collate_path)) {
        # TODO: delete this if this does nothing!
        colData.Rds <- list()
        colData.Rds$df.anno <- df.anno
        colData.Rds$df.files <- df.files
        colData.Rds$bam_path <- bam_path
        colData.Rds$sw_path <- sw_path        
    }
}

.infobox_update_se <- function(se, path) {
    ui_infobox_expr(ifelse(
        is(se, "NxtSE"), 2, ifelse(
            is_valid(path) && file.exists(file.path(path,"colData.Rds")),
            1,0)))
}

.server_expr_load_alert_success <- function(session, collate_path) {
    sendSweetAlert(
        session = session,
        title = paste("Experiment Loaded successfully from", 
            collate_path),
        type = "success"
    )
}

.server_expr_load_alert_fail <- function(session, collate_path) {
    sendSweetAlert(
        session = session,
        title = paste("No valid experiment found at", 
            collate_path),
        type = "error"
    )
}

.server_expr_ref_load_success <- function(session, ref_path) {
    sendSweetAlert(
        session = session,
        title = paste("Reference loaded successfully from", 
            ref_path),
        type = "success"
    )
}

.server_expr_ref_load_fail <- function(session, ref_path) {
    sendSweetAlert(
        session = session,
        title = paste("Reference loading failed from", 
            ref_path),
        type = "error"
    )
}

.makeSE_sweetalert_finish <- function(session) {
    sendSweetAlert(
        session = session,
        title = "NxtSE object loaded successfully",
        type = "success"
    )
}

.makeSE_sweetalert_error <- function(session) {
    sendSweetAlert(
        session = session,
        title = "Error encountered loading NxtSE object",
        type = "error"
    )
}

.load_NxtSE_sweetalert_finish <- function(session) {
    sendSweetAlert(
        session = session,
        title = "Successfully loaded NxtSE from RDS",
        type = "success"
    )
}

.load_NxtSE_sweetalert_error <- function(session) {
    sendSweetAlert(
        session = session,
        title = "Error encountered loading NxtSE from RDS",
        type = "error"
    )
}