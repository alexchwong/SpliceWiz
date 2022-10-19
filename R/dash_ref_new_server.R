# Module for new SpliceWiz reference
server_ref_new <- function(id, refresh_tab, volumes, get_memmode_reactive) {
    moduleServer(id, function(input, output, session) {
        # Instantiate settings
        settings_newref <- setreactive_newref()
        
        # Choose path to store reference
        observe({  
            shinyDirChoose(input, "dir_reference_path", 
                roots = volumes(), session = session)
            # display reference path
            output$txt_reference_path <- renderText({
                validate(need(input$dir_reference_path, 
                    "Please select reference path"))
                settings_newref$newref_path <- 
                    parseDirPath(volumes(), input$dir_reference_path)
            })
        })
        
        # Choose genome FASTA
        observe({
            shinyFileChoose(input, "file_genome", 
                roots = volumes(), 
                session = session, filetypes = c("fa", "fasta", "gz"))
            if(!is.null(input$file_genome)){
                file_selected<-parseFilePaths(volumes(), input$file_genome)
                settings_newref$newref_fasta <- 
                    as.character(file_selected$datapath)
                output$txt_genome <- renderText(
                    as.character(file_selected$datapath))
            }
        })
        
        # Choose annotation GTF
        observe({  
            shinyFileChoose(input, "file_gtf", roots = volumes(), 
                session = session, filetypes = c("gtf", "gz"))
            if(!is.null(input$file_gtf)){
                file_selected<-parseFilePaths(volumes(), input$file_gtf)
                settings_newref$newref_gtf <- 
                    as.character(file_selected$datapath)
                output$txt_gtf <- renderText(
                    as.character(file_selected$datapath))
            }
        })

        # Choose mappability file
        observe({  
            shinyFileChoose(input, "file_mappa", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz", "rds"))
            if(!is.null(input$file_mappa)){
                # updateSelectInput(session = session, 
                    # inputId = "newref_genome_type", 
                    # choices = c("(custom)", "hg38", "mm10", "hg19", "mm9"),
                    # selected = "(custom)"
                # )
                file_selected <- parseFilePaths(volumes(), input$file_mappa)
                settings_newref$newref_mappa <- 
                    as.character(file_selected$datapath)
            }
        })
        observeEvent(settings_newref$newref_mappa, {
            output$txt_mappa <- renderText(settings_newref$newref_mappa)
        })
        observeEvent(input$clear_mappa, {
            req(input$clear_mappa)
            settings_newref$newref_mappa <- ""
        })    

        # Choose non-polyA BED file
        observe({  
            shinyFileChoose(input, "file_NPA", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz"))
            if(!is.null(input$file_NPA)){
                # updateSelectInput(session = session, 
                    # inputId = "newref_genome_type", 
                    # choices = c("(custom)", "hg38", "mm10", "hg19", "mm9"),
                    # selected = "(custom)"
                # )
                file_selected <- parseFilePaths(volumes(), input$file_NPA)
                settings_newref$newref_NPA <- 
                    as.character(file_selected$datapath)
            }
        })
        observeEvent(settings_newref$newref_NPA, {
            output$txt_NPA <- renderText(settings_newref$newref_NPA)    
        })
        observeEvent(input$clear_NPA, {
            req(input$clear_NPA)
            settings_newref$newref_NPA <- ""
        })
        
        # Choose Blacklist BED file
        observe({  
            shinyFileChoose(input, "file_bl", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz"))
            if(!is.null(input$file_bl)){
                file_selected <- parseFilePaths(volumes(), input$file_bl)
                settings_newref$newref_bl <- 
                    as.character(file_selected$datapath)
            }
        })
        observeEvent(settings_newref$newref_bl, {
            output$txt_bl <- renderText(settings_newref$newref_bl)
        })
        observeEvent(input$clear_bl, {
            req(input$clear_bl)
            settings_newref$newref_bl <- ""
        })
        
        # Choose genome_type
        observeEvent(input$newref_genome_type, {
            req(input$newref_genome_type)

            gt <- input$newref_genome_type
            valid_gt_options <- c("hg38", "hg19", "mm10", "mm9")
            withProgress(message = "Retrieving Mappability resource", value = 0,
            {
                if(gt %in% valid_gt_options) {
                    settings_newref$newref_NPA <- getNonPolyARef(gt)
                    settings_newref$newref_mappa <- 
                        .sw_dash_get_mappa(gt)
                } else if(gt == "(custom)") {
            # do nothing. This allows user to first select the default 
            #   and then change to user-defined files
                } else {
                    settings_newref$newref_NPA <- ""
                    settings_newref$newref_mappa <- ""
                }
            })

            settings_newref$ui_newref_genome_type <- input$newref_genome_type
        })

        # This block runs when the tab is refreshed
        #   or if `release` or `species` is changed
        observeEvent({
            list(
                refresh_tab(), input$release, input$species
            )
        }, {
            req(refresh_tab())
            msg <- 'Fetching from Ensembl FTP'
            if(is_valid(input$release) & is_valid(input$species)) {
                withProgress(message = msg, value = 0, {
                    updateSelectInput(session = session, 
                        inputId = "fasta", 
                        choices = c("", .refresh_genome(
                            input$release, input$species
                        ))
                    )
                    updateSelectInput(session = session, 
                        inputId = "gtf", 
                        choices = c("", .refresh_gtf(
                            input$release, input$species
                        ))
                    )
                })
            } else if(is_valid(input$release)) {
                withProgress(message = msg, value = 0, {
                    updateSelectInput(session = session, 
                        inputId = "species", 
                        choices = c("", .refresh_species(
                            input$release
                        ))
                    )
                })            
            } else {
                withProgress(message = msg, value = 0, {
                    updateSelectInput(session = session, 
                        inputId = "release", 
                        choices = c("", .refresh_releases())
                    )
                })
            }
        })
        
        # Choose FASTA from Ensembl
        observeEvent(input$fasta, {
            req(input$fasta)
            settings_newref$newref_fasta <- paste0(
                # "https://ftp.ensembl.org/pub/",
                "ftp://ftp.ensembl.org/pub/",
                "release-", as.character(isolate(input$release)),
                "/fasta/", isolate(input$species), "/dna/",
                input$fasta
            )
            output$txt_genome <- renderText(settings_newref$newref_fasta)
        })

        # Choose GTF from Ensembl
        observeEvent(input$gtf, {
            req(input$gtf)
            settings_newref$newref_gtf <- paste0(
                # "https://ftp.ensembl.org/pub/",
                "ftp://ftp.ensembl.org/pub/",
                "release-", as.character(isolate(input$release)),
                "/gtf/", isolate(input$species), "/",
                input$gtf
            )
            output$txt_gtf <- renderText(settings_newref$newref_gtf)
        })
        
        # Click Build Reference
        observeEvent(input$buildRef, {
            args <- list(
                reference_path  = settings_newref$newref_path, 
                fasta           = settings_newref$newref_fasta, 
                gtf             = settings_newref$newref_gtf,
                genome_type     = input$newref_genome_type, 
                nonPolyARef     = settings_newref$newref_NPA, 
                MappabilityRef  = settings_newref$newref_mappa,
                BlacklistRef    = settings_newref$newref_bl
            )

            args <- Filter(is_valid, args)
            if(!("reference_path" %in% names(args))) {
                output$refStatus <- renderText({ "Reference path not set" })
            } else if(!any(c("fasta") %in% names(args))) {
                output$refStatus <- renderText({ "Genome not provided" })        
            } else if(!any(c("gtf") %in% names(args))) {
                output$refStatus <- renderText("Gene annotations not provided")
            } else {        
                # Copy MappabilityRef into target directory
                if(
                        "MappabilityRef" %in% names(args) && 
                        file.exists(args$MappabilityRef)
                ) {
                    mappa_base <- basename(args$MappabilityRef)
                    new_mappa_path <- file.path(args$reference_path, 
                        "Mappability")
                    new_mappa_file <- file.path(new_mappa_path, mappa_base)
                    
                    if(!dir.exists(new_mappa_path)) dir.create(new_mappa_path)
                    
                    if(dir.exists(new_mappa_path))
                        file.copy(args$MappabilityRef, new_mappa_file)
                    
                    if(file.exists(new_mappa_file)) 
                        args$MappabilityRef <- new_mappa_file
                }
                args$lowMemoryMode <- get_memmode_reactive()
                withProgress(message = 'Building Reference', value = 0, {
                    do.call(buildRef, args)
                })
                # If successfully created, load this reference automatically
                if(file.exists(
                    file.path(settings_newref$newref_path, "settings.Rds")
                )) {
                    sendSweetAlert(
                        session = session,
                        title = "Reference Build complete!",
                        type = "success"
                    )           
                } else {
                    sendSweetAlert(
                        session = session,
                        title = paste("Reference Build failed.",
                            "An error must have occurred"),
                        type = "error"
                    )               
                }
            }
        })
            
        # clearNewRef Button
        observeEvent(input$clearNewRef, {
            settings_newref <- setreactive_newref()
            output$txt_reference_path <- 
                renderText("Please select reference path")
            output$txt_genome <- renderText("")
            output$txt_gtf <- renderText("")
            output$txt_mappa <- renderText("")
            output$txt_NPA <- renderText("")
            output$txt_bl <- renderText("")
            updateSelectInput(session = session, 
                inputId = "fasta", 
                choices = c("")
            )
            updateSelectInput(session = session, 
                inputId = "gtf", 
                choices = c("")
            )
            updateSelectInput(session = session, 
                inputId = "species", 
                choices = c("")
            )
            updateSelectInput(session = session, 
                inputId = "release", 
                choices = c("")
            )
        })
        
        observeEvent(input$load_ref_example, {
            if(!dir.exists(file.path(tempdir(), "Reference")))
                dir.create(file.path(tempdir(), "Reference"))
            output$txt_reference_path <- renderText({
                settings_newref$newref_path <- file.path(tempdir(), "Reference")
            })
            settings_newref$newref_fasta <- NxtIRFdata::chrZ_genome()
            settings_newref$newref_gtf <- NxtIRFdata::chrZ_gtf()
            
            output$txt_genome <- renderText(settings_newref$newref_fasta)
            output$txt_gtf <- renderText(settings_newref$newref_gtf)
        })
        
        return(settings_newref)
    })
}

.refresh_releases <- function() {
    test <- XML::getHTMLLinks("http://ftp.ensembl.org/pub")
    test <- test[grepl("release-", test)]
    test <- test[grepl("/", test)]
    int_release <- tstrsplit(test, split="-")[[2]]
    int_release <- as.integer(sub("/","",int_release))
    int_release <- sort(int_release, decreasing = TRUE)
    return(int_release[int_release > 46])
}

.refresh_species <- function(release) {
    test_genome <- XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/fasta/"
    ))
    test_gtf <- XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/gtf/"
    ))
    species <- union(test_genome, test_gtf)
    species <- species[!grepl("..", species, fixed = TRUE)]
    species <- sub("/","",species)
    if(all(c("homo_sapiens", "mus_musculus") %in% species)) {
        species <- unique(c(
            "homo_sapiens",
            "mus_musculus",
            sort(species)
            )
        )        
    } else {
        species <- sort(species)
    }
    return(species)
}

.refresh_genome <- function(release, species) {
    test_genome <- XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/fasta/",
        species, "/dna/"           
    ))
    test_genome <- test_genome[
        grepl("toplevel", test_genome, fixed = TRUE) |
        grepl("primary", test_genome, fixed = TRUE)
    ]
    test_genome
}

.refresh_gtf <- function(release, species) {
    test_gtf <- XML::getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/gtf/",
        species        
    ))
    test_gtf <- test_gtf[
        grepl("gtf.gz", test_gtf, fixed = TRUE) &
        !grepl("abinitio", test_gtf, fixed = TRUE) &
        !grepl("scaff", test_gtf, fixed = TRUE) &
        !grepl(".chr.", test_gtf, fixed = TRUE)
    ]
    test_gtf
}

.sw_dash_get_mappa <- function(genome_type, path = tempdir()) {
    temp <- ""
    tryCatch({
        temp <- get_mappability_exclusion(
            genome_type, as_type = "bed", path)
    }, error = function(e) {
        message(e)
        temp <- ""
    })
    if(file.exists(temp)) {
        return(temp)
    } else {
        return("")
    }
}