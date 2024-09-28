# Module for new SpliceWiz reference
server_ref_new <- function(id, refresh_tab, volumes, get_memmode_reactive) {
    moduleServer(id, function(input, output, session) {
        # Instantiate settings
        settings_newref <- setreactive_newref()
        
        # Choose path to store reference
        observe({  
            shinyDirChoose(input, "dir_reference_path", 
                roots = volumes(), session = session)
            settings_newref$newref_path <- 
                parseDirPath(volumes(), input$dir_reference_path)
        })
        
        # Choose genome FASTA
        observe({
            shinyFileChoose(input, "file_genome", 
                roots = volumes(), 
                session = session, filetypes = c("fa", "fasta", "gz"))
            req(input$file_genome)
            file_selected <- parseFilePaths(volumes(), input$file_genome)
            settings_newref$newref_fasta <- 
                as.character(file_selected$datapath)
        })
        
        # Choose annotation GTF
        observe({  
            shinyFileChoose(input, "file_gtf", roots = volumes(), 
                session = session, filetypes = c("gtf", "gz"))
            req(input$file_gtf)
            file_selected<-parseFilePaths(volumes(), input$file_gtf)
            settings_newref$newref_gtf <- 
                as.character(file_selected$datapath)
        })

        # Choose mappability file
        observe({  
            shinyFileChoose(input, "file_mappa", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz", "rds"))
            req(input$file_mappa)
            file_selected <- parseFilePaths(volumes(), input$file_mappa)
            settings_newref$newref_mappa <- 
                as.character(file_selected$datapath)
        })

        observeEvent(input$clear_mappa, {
            req(input$clear_mappa)
            settings_newref$newref_mappa <- ""
        })    

        # Choose non-polyA BED file
        observe({  
            shinyFileChoose(input, "file_NPA", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz"))
            req(input$file_NPA)
            file_selected <- parseFilePaths(volumes(), input$file_NPA)
            settings_newref$newref_NPA <- 
                as.character(file_selected$datapath)
        })

        observeEvent(input$clear_NPA, {
            req(input$clear_NPA)
            settings_newref$newref_NPA <- ""
        })
        
        # Choose Blacklist BED file
        observe({  
            shinyFileChoose(input, "file_bl", roots = volumes(), 
                session = session, filetypes = c("bed", "txt", "gz"))
            req(input$file_bl)
            file_selected <- parseFilePaths(volumes(), input$file_bl)
            settings_newref$newref_bl <- 
                as.character(file_selected$datapath)
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
            if(gt %in% valid_gt_options) {
                withProgress(
                        message = "Retrieving Mappability resource", value = 0,
                {
                    settings_newref$newref_NPA <- getNonPolyARef(gt)
                    settings_newref$newref_mappa <- .sw_dash_get_mappa(gt)
                })
            }

            if(gt %in% c("hg38", "hg19")) {
                updateSelectInput(session = session, 
                    inputId = "newref_species_GO", 
                    choices = c("(none)", "Homo sapiens", "Mus musculus"),
                    selected = "Homo sapiens"
                )            
            } else if(gt %in% c("mm10", "mm9")) {
                updateSelectInput(session = session, 
                    inputId = "newref_species_GO", 
                    choices = c("(none)", "Homo sapiens", "Mus musculus"),
                    selected = "Mus musculus"
                )   
            }
        })

        # Refresh release
        observeEvent(refresh_tab(), {
            # Check Ensembl is working
            if(!is_valid(settings_newref$ensemblWorking)) {
                    withProgress(
                            message = "Testing Ensembl FTP is working", value = 0, 
                    {
                        test <- .check_ensembl()
                    })
                if(is(test, "list")) {
                    settings_newref$ensemblWorking <- test                
                } else {
                    withProgress(
                            message = "Retrieving data from Ensembl FTP", value = 0, 
                    {
                        settings_newref$availRelease <- c("", .refresh_releases())
                    })                                    
                }
            }
        })
        
        # Refresh species
        observeEvent(input$release, {
            req(is_valid(input$release))
            withProgress(
                    message = "Retrieving data from Ensembl FTP", value = 0, 
            {
                settings_newref$availSpecies <- c("", 
                    .refresh_species(input$release))
            })        
        })

        # Refresh genome / GTF
        observeEvent(input$species, {
            req(is_valid(input$release))
            req(is_valid(input$species))
            req(input$species %in% settings_newref$availSpecies)
            withProgress(
                    message = "Retrieving data from Ensembl FTP", value = 0, 
            {
                settings_newref$availGenome <- c("", .refresh_genome(
                    input$release, input$species
                ))
                settings_newref$availGTF <- c("", .refresh_gtf(
                    input$release, input$species
                ))
            })
        })
        
        # Choose FASTA from Ensembl
        observeEvent(input$fasta, {
            req(is_valid(input$release))
            req(is_valid(input$species))
            req(input$species %in% settings_newref$availSpecies)
            req(is_valid(input$fasta))
            settings_newref$newref_fasta <- paste0(
                # "https://ftp.ensembl.org/pub/",
                "ftp://ftp.ensembl.org/pub/",
                "release-", as.character(isolate(input$release)),
                "/fasta/", isolate(input$species), "/dna/",
                input$fasta
            )
        })

        # Choose GTF from Ensembl
        observeEvent(input$gtf, {
            req(is_valid(input$release))
            req(is_valid(input$species))
            req(input$species %in% settings_newref$availSpecies)
            req(is_valid(input$gtf))
            settings_newref$newref_gtf <- paste0(
                # "https://ftp.ensembl.org/pub/",
                "ftp://ftp.ensembl.org/pub/",
                "release-", as.character(isolate(input$release)),
                "/gtf/", isolate(input$species), "/",
                input$gtf
            )
        })

        observeEvent(input$newref_species_GO, {
            req(input$newref_species_GO)
            settings_newref$newref_speciesGO <- input$newref_species_GO
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
                BlacklistRef    = settings_newref$newref_bl,
                ontologySpecies = settings_newref$newref_speciesGO
            )

            args <- Filter(is_valid, args)
            if(!("reference_path" %in% names(args))) {
                output$refStatus <- renderText({ "Reference path not set" })
            } else if(!any(c("fasta") %in% names(args))) {
                output$refStatus <- renderText({ "Genome not provided" })        
            } else if(!any(c("gtf") %in% names(args))) {
                output$refStatus <- renderText("Gene annotations not provided")
            } else {        
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
            settings_newref$newref_path <- ""
            settings_newref$newref_fasta <- ""
            settings_newref$newref_gtf <- ""
            settings_newref$newref_NPA <- ""
            settings_newref$newref_mappa <- ""
            settings_newref$newref_bl <- ""

            settings_newref$availSpecies <- ""
            settings_newref$availGenome <- ""
            settings_newref$availGTF <- ""

            updateSelectInput(session = session, 
                inputId = "newref_genome_type", 
                choices = c("(none)", "hg38", "mm10", "hg19", "mm9"),
                selected = "(none)"
            )

            updateSelectInput(session = session, 
                inputId = "newref_species_GO", 
                choices = c("(none)", "Homo sapiens", "Mus musculus"),
                selected = "(none)"
            )
        })
        
        observeEvent(input$load_ref_example, {
            if(!dir.exists(file.path(tempdir(), "Reference")))
                dir.create(file.path(tempdir(), "Reference"))
            updateSelectInput(session = session, 
                inputId = "newref_species_GO", 
                choices = c("(none)", "Homo sapiens", "Mus musculus"),
                selected = "Homo sapiens"
            )
            settings_newref$newref_path <- file.path(tempdir(), "Reference")
            settings_newref$newref_fasta <- NxtIRFdata::chrZ_genome()
            settings_newref$newref_gtf <- NxtIRFdata::chrZ_gtf()
            output$txt_demo <- renderText({
                paste(
                    "Demo reference settings loaded.",
                    "Click `Build Reference` to build demo SpliceWiz reference"
                )
            })
        })

    # Output displays
        output$txt_reference_path <- renderText({
            validate(need(settings_newref$newref_path, 
                "Please select reference path"))
            settings_newref$newref_path
        })
        output$txt_genome <- renderText(settings_newref$newref_fasta)
        output$txt_gtf <- renderText(settings_newref$newref_gtf)
        output$txt_mappa <- renderText(settings_newref$newref_mappa)
        output$txt_NPA <- renderText(settings_newref$newref_NPA)    
        output$txt_bl <- renderText(settings_newref$newref_bl)
    
    # Drop-down options based on settings
        observeEvent(settings_newref$ensemblWorking, {
            if(is(settings_newref$ensemblWorking, "list")) {
                errMsg <- settings_newref$ensemblWorking$error
                sendSweetAlert(
                    session = session,
                    title = paste("Unable to access Ensembl -", errMsg),
                    type = "error"
                )   
            }
        })
        observeEvent(settings_newref$availRelease, {
            updateSelectInput(session = session, 
                inputId = "release", 
                choices = c("", settings_newref$availRelease)
            )        
        })
        observeEvent(settings_newref$availSpecies, {
            updateSelectInput(session = session, 
                inputId = "species", 
                choices = c("", settings_newref$availSpecies)
            )        
        })
        observeEvent(settings_newref$availGenome, {
            updateSelectInput(session = session, 
                inputId = "fasta", 
                choices = c("", settings_newref$availGenome)
            )        
        })
        observeEvent(settings_newref$availGTF, {
            updateSelectInput(session = session, 
                inputId = "gtf", 
                choices = c("", settings_newref$availGTF)
            )        
        })
        
        return(settings_newref)
    })
}

.check_ensembl <- function() {
    url <- "http://ftp.ensembl.org/pub"
    out <- tryCatch(
        {
            GET(url, timeout(10)) %>% read_html() %>% html_nodes("a") %>%html_attr("href")
        },
        error=function(cond) {
            errList <- list(error = cond)
            return(errList)
        }
    )
    return(out)
}

.refresh_releases <- function() {
    test <- getHTMLLinks("http://ftp.ensembl.org/pub")
    test <- test[grepl("release-", test)]
    test <- test[grepl("/", test)]
    int_release <- tstrsplit(test, split="-")[[2]]
    int_release <- as.integer(sub("/","",int_release))
    int_release <- sort(int_release, decreasing = TRUE)
    return(int_release[int_release > 46])
}

.refresh_species <- function(release) {
    test_genome <- getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/fasta/"
    ))
    test_genome <- test_genome[grepl("/", test_genome)]
    test_gtf <- getHTMLLinks(paste0(
        "http://ftp.ensembl.org/pub/",
        "release-",
        as.character(release),
        "/gtf/"
    ))
    test_gtf <- test_gtf[grepl("/", test_gtf)]
    species <- union(test_genome, test_gtf)
    species <- species[!grepl("..", species, fixed = TRUE)]
    species <- species[!grepl("release", species, fixed = TRUE)]
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
    test_genome <- getHTMLLinks(paste0(
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
    test_gtf <- getHTMLLinks(paste0(
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