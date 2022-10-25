server_DE <- function(
        id, refresh_tab, volumes, get_threads,
        get_se, get_filters, get_selected_diag, get_selected_volc
) {
    moduleServer(id, function(input, output, session) {
        settings_DE <- setreactive_DE()

        observeEvent(refresh_tab(), {
            # output$warning_DE = renderText({
                # validate(need(get_se(), 
                    # "Please load experiment via 'Experiment' tab"))
                
                # "Experiment Loaded"
            # })
            req(get_se())
            colcat <- .get_discrete_cats(get_se())
            updateSelectInput(session = session, inputId = "variable_DE", 
                choices = c("(none)", colcat), selected = "(none)")
            updateSelectInput(session = session, inputId = "batch1_DE", 
                choices = c("(none)", colcat), selected = "(none)")
            updateSelectInput(session = session, inputId = "batch2_DE", 
                choices = c("(none)", colcat), selected = "(none)")
        })
        
        observeEvent(input$DT_DE_rows_all, {
            settings_DE$DT_DE_rows_all <- input$DT_DE_rows_all
        })
        
        observeEvent(input$DT_DE_rows_selected, {
            settings_DE$DT_DE_rows_selected <- input$DT_DE_rows_selected
        })
        
        # Filter import
        observeEvent(get_filters(), {
            settings_DE$filters <- get_filters()
        })
        
        observeEvent(input$variable_DE, {
            req(input$variable_DE)
            req(input$variable_DE != "(none)")

            colData <- colData(get_se())
            
            if(is(colData[,input$variable_DE], "factor")) {
                col_levels <- levels(colData[,input$variable_DE])
                updateSelectInput(session = session, inputId = "nom_DE", 
                    choices = c("(none)", col_levels), 
                    selected = "(none)")
                updateSelectInput(session = session, inputId = "denom_DE", 
                    choices = c("(none)", col_levels),
                    selected = "(none)")
                if(
                        is_valid(settings_DE$nom_DE) && 
                        settings_DE$nom_DE %in% col_levels
                ) {
                    updateSelectInput(session = session, inputId = "nom_DE", 
                        selected = settings_DE$nom_DE)
                }
                if(
                        is_valid(settings_DE$nom_DE) && 
                        settings_DE$denom_DE %in% col_levels
                ) {
                    updateSelectInput(session = session, inputId = "denom_DE", 
                        selected = settings_DE$denom_DE)
                }
            } else {
                updateSelectInput(session = session, inputId = "nom_DE", 
                    choices = c("(time series)"), 
                    selected = "(time series)")
                updateSelectInput(session = session, inputId = "denom_DE", 
                    choices = c("(time series denom)"),
                    selected = "(time series denom)")
            }
        })

        observeEvent(settings_DE$method, {
            req(settings_DE$method)
            updateSelectInput(session = session, inputId = "method_DE", 
                selected = settings_DE$method)
        })
        observeEvent(settings_DE$DE_Var, {
            req(settings_DE$DE_Var)
            updateSelectInput(session = session, inputId = "variable_DE", 
                selected = settings_DE$DE_Var)
        })
        observeEvent(settings_DE$nom_DE, {
            req(settings_DE$nom_DE)
            updateSelectInput(session = session, inputId = "nom_DE_DE", 
                selected = settings_DE$nom_DE)
        })
        observeEvent(settings_DE$denom_DE, {
            req(settings_DE$denom_DE)
            updateSelectInput(session = session, inputId = "denom_DE", 
                selected = settings_DE$denom_DE)
        })
        observeEvent(settings_DE$batchVar1, {
            req(settings_DE$batchVar1)
            updateSelectInput(session = session, inputId = "batch1_DE", 
                selected = settings_DE$batchVar1)
        })
        observeEvent(settings_DE$batchVar2, {
            req(settings_DE$batchVar2)
            updateSelectInput(session = session, inputId = "batch2_DE", 
                selected = settings_DE$batchVar2)
        })
        
        observeEvent(settings_DE$IRmode_DE, {
            req(settings_DE$IRmode_DE)
            shinyWidgets::updateRadioGroupButtons(
                session = session, inputId = "de_IRmode", 
                selected = settings_DE$IRmode_DE)
        })
        
        observe({
            output$warning_DE <- renderText({
                validate(need(get_se(), 
                    "Please load experiment via 'Experiment' tab"))
                validate(need(is_valid(input$variable_DE),
                    "Variable for DE needs to be defined"))
                validate(need(is_valid(input$nom_DE), 
                    "Nominator for DE Variable needs to be defined"))
                validate(need(is_valid(input$denom_DE), 
                    "Denominator for DE Variable needs to be defined"))
                validate(need(input$nom_DE != "(time series)" ||
                        input$method_DE == "DESeq2", 
                    "Time series analysis can only be performed using DESeq2"))
                validate(need(input$denom_DE != input$nom_DE, 
                    "Denominator and Nominator must be different"))
                "Ready to run differential analysis"
            })
        })

        observeEvent(input$perform_DE, {
            req(get_se())
            output$warning_DE <- renderText({
                validate(need(is_valid(input$variable_DE),
                    "Variable for DE needs to be defined"))
                validate(need(is_valid(input$nom_DE), 
                    "Nominator for DE Variable needs to be defined"))
                validate(need(is_valid(input$denom_DE), 
                    "Denominator for DE Variable needs to be defined"))
                validate(need(input$nom_DE != "(time series)" ||
                        input$method_DE == "DESeq2", 
                    "Time series analysis can only be performed using DESeq2"))
                validate(need(input$denom_DE != input$nom_DE, 
                    "Denominator and Nominator must be different"))
                "Running differential analysis"
            })
            req(is_valid(input$variable_DE))
            req(is_valid(input$nom_DE))
            req(is_valid(input$denom_DE))
            req(input$denom_DE != input$nom_DE)
            req(input$nom_DE != "(time series)" ||
                input$method_DE == "DESeq2")

            rowData <- as.data.frame(rowData(get_se()))
            colData <- as.data.frame(colData(get_se()))

            settings_DE$DE_Var <- input$variable_DE
            settings_DE$nom_DE <- input$nom_DE
            settings_DE$denom_DE <- input$denom_DE

            if(
                    input$batch1_DE != "(none)" & 
                    input$batch1_DE != input$variable_DE
            ) {
                settings_DE$batchVar1 <- input$batch1_DE
            } else {
                settings_DE$batchVar1 <- ""
                updateSelectInput(session = session, inputId = "batch1_DE", 
                    selected = "(none)")
            }
            if(
                    input$batch2_DE != "(none)" & 
                    input$batch2_DE != input$variable_DE & 
                    input$batch2_DE != input$batch1_DE
            ) {
                settings_DE$batchVar2 <- input$batch2_DE
            } else {
                settings_DE$batchVar2 <- ""
                updateSelectInput(session = session, inputId = "batch2_DE", 
                    selected = "(none)")
            }
            settings_DE$IRmode_DE <- input$de_IRmode
            
            req(input$method_DE)
            settings_DE$method <- input$method_DE

            if(settings_DE$method == "DESeq2") {
                withProgress(message = 'Running DESeq2...', value = 0, {
                    if(settings_DE$nom_DE != "(time series)") {
                        res.ASE <- ASE_DESeq(
                            se = get_se(), 
                            test_factor = settings_DE$DE_Var, 
                            test_nom = settings_DE$nom_DE, 
                            test_denom = settings_DE$denom_DE,
                            batch1 = settings_DE$batchVar1, 
                            batch2 = settings_DE$batchVar2,
                            n_threads = get_threads(),
                            IRmode = settings_DE$IRmode_DE
                        )                    
                    } else {
                        res.ASE <- ASE_DESeq(
                            se = get_se(), 
                            test_factor = settings_DE$DE_Var, 
                            batch1 = settings_DE$batchVar1, 
                            batch2 = settings_DE$batchVar2,
                            n_threads = get_threads(),
                            IRmode = settings_DE$IRmode_DE
                        )
                    }
                    
                    if(!input$adjP_DE) {
                        setorderv(res.ASE, "pvalue")
                    } else {
                        setorderv(res.ASE, "padj")            
                    }
                    settings_DE$res <- as.data.frame(res.ASE)
                })
            } else if(settings_DE$method == "limma") {
                withProgress(message = 'Running limma...', value = 0, {
                    res.ASE <- ASE_limma(
                        se = get_se(), 
                        test_factor = settings_DE$DE_Var, 
                        test_nom = settings_DE$nom_DE, 
                        test_denom = settings_DE$denom_DE,
                        batch1 = settings_DE$batchVar1, 
                        batch2 = settings_DE$batchVar2,
                        IRmode = settings_DE$IRmode_DE
                    )
                    if(!input$adjP_DE) {
                        setorderv(res.ASE, "P.Value")
                    } else {
                        setorderv(res.ASE, "adj.P.Val")
                    }
                })
            } else if(settings_DE$method == "DoubleExpSeq") {
                withProgress(message = 'Running DoubleExpSeq...', value = 0, {
                    res.ASE <- ASE_DoubleExpSeq(
                        se = get_se(), 
                        test_factor = settings_DE$DE_Var, 
                        test_nom = settings_DE$nom_DE, 
                        test_denom = settings_DE$denom_DE,
                        IRmode = settings_DE$IRmode_DE
                    )
                    if(!input$adjP_DE) {
                        setorderv(res.ASE, "P.Value")
                    } else {
                        setorderv(res.ASE, "adj.P.Val")
                    }
                })
            } else if(settings_DE$method == "satuRn") {
                withProgress(message = 'Running satuRn...', value = 0, {
                    res.ASE <- ASE_satuRn(
                        se = get_se(), 
                        test_factor = settings_DE$DE_Var, 
                        test_nom = settings_DE$nom_DE, 
                        test_denom = settings_DE$denom_DE,
                        batch1 = settings_DE$batchVar1, 
                        batch2 = settings_DE$batchVar2,
                        IRmode = settings_DE$IRmode_DE,
                        n_threads = get_threads()
                    )
                    if(!input$adjP_DE) {
                        setorderv(res.ASE, "pval")
                    } else {
                        setorderv(res.ASE, "regular_FDR")
                    }
                })
            }
            
            # Allow filtering by EventType: factorize it
            if(
                    is(res.ASE$EventType, "character") && 
                    length(unique(res.ASE$EventType)) > 1
            ) {
                res.ASE$EventType <- factor(res.ASE$EventType)
            }
            
            settings_DE$res <- as.data.frame(res.ASE)
            output$warning_DE <- renderText({"Finished"})

            req(settings_DE$res)
            # save settings for current res
            settings_DE$res_settings$method <- settings_DE$method
            settings_DE$res_settings$DE_Var <- settings_DE$DE_Var
            settings_DE$res_settings$nom_DE <- settings_DE$nom_DE
            settings_DE$res_settings$denom_DE <- settings_DE$denom_DE
            settings_DE$res_settings$batchVar1 <- settings_DE$batchVar1
            settings_DE$res_settings$batchVar2 <- settings_DE$batchVar2
            settings_DE$res_settings$IRmode_DE <- settings_DE$IRmode_DE
            settings_DE$res_settings$BuildVersion <- ASE_version
        })

        observeEvent(settings_DE$res, {
            req(settings_DE$res)
            output$DT_DE <- DT::renderDataTable(
                DT::datatable(
                    settings_DE$res,
                    class = 'cell-border stripe',
                    rownames = FALSE,
                    filter = 'top'
                )
            )
        })
        
        observe({
            shinyFileSave(input, "save_DE", roots = volumes(), 
                session = session, filetypes = c("Rds"))
        })
        observeEvent(input$save_DE, {
            req(settings_DE$res)
            req(length(settings_DE$res_settings) > 0)
            req("BuildVersion" %in% names(settings_DE$res_settings))

            selectedfile <- parseSavePath(volumes(), input$save_DE)
            req(selectedfile$datapath)

            save_DE <- list(
                res = settings_DE$res, 
                settings = settings_DE$res_settings, 
                filters = settings_DE$filters)
                
            output$warning_DE <- renderText({
                validate(need(
                    tryCatch({
                        saveRDS(save_DE,selectedfile$datapath)
                        TRUE
                        }, error = function(e) {                            
                            message(e)
                            FALSE
                        }
                    ),
                    "An error occurred attempting to save to file."
                ))
            })
        })
        
        observe({
            shinyFileChoose(input, "load_DE", roots = volumes(), 
                session = session, filetype = c("Rds"))
        })
        observeEvent(input$load_DE, {
            req(input$load_DE)
            file_selected <- parseFilePaths(volumes(), input$load_DE)
            req(file_selected$datapath)

            load_DE <- NULL
            retVal <- tryCatch({
                load_DE <- readRDS(as.character(file_selected$datapath))
                TRUE
            }, error = function(e) {                            
                message(e)
                FALSE
            })

            output$warning_DE <- renderText({
                validate(need(retVal,
                    "An error occurred attempting to load from file."
                ))
            })
            
            req(all(c("res", "settings") %in% names(load_DE)))
            # check all parameters exist in colData(se)
            output$warning_DE <- renderText({
                validate(need(get_se(), 
                    "Please load experiment via 'Experiment' tab"))
                validate(need(
                        "BuildVersion" %in% names(load_DE$settings) &&
                        load_DE$settings$BuildVersion >= ASE_version, 
                    paste(
                        "This differential expression Rds was produced with",
                        "SpliceWiz version 0.99.4 or below."
                    )))
                
                "Experiment Loaded"
            })
            req(get_se())
            colData <- colData(get_se())
            req(load_DE$settings$DE_Var %in% colnames(colData))
            req(!is_valid(load_DE$settings$batchVar1) || 
                load_DE$settings$batchVar1 %in% colnames(colData))
            req(!is_valid(load_DE$settings$batchVar2) || 
                load_DE$settings$batchVar2 %in% colnames(colData))
            req(any(unlist(colData[,load_DE$settings$DE_Var]) == 
                load_DE$settings$nom_DE))
            req(any(unlist(colData[,load_DE$settings$DE_Var]) == 
                load_DE$settings$denom_DE))
            req(load_DE$settings$method %in% 
                c("DESeq2", "limma", "DoubleExpSeq", "satuRn"))
            req(load_DE$settings$IRmode_DE %in% 
                c("all", "annotated", "annotated_binary"))

            settings_DE$res <- load_DE$res
            settings_DE$res_settings$method <- load_DE$settings$method
            settings_DE$res_settings$DE_Var <- load_DE$settings$DE_Var
            settings_DE$res_settings$nom_DE <- load_DE$settings$nom_DE
            settings_DE$res_settings$denom_DE <- load_DE$settings$denom_DE
            settings_DE$res_settings$batchVar1 <- load_DE$settings$batchVar1
            settings_DE$res_settings$batchVar2 <- load_DE$settings$batchVar2
            settings_DE$res_settings$IRmode_DE <- load_DE$settings$IRmode_DE
            settings_DE$res_settings$BuildVersion <- 
                load_DE$settings$BuildVersion

            settings_DE$method <- settings_DE$res_settings$method
            settings_DE$DE_Var <- settings_DE$res_settings$DE_Var
            settings_DE$nom_DE <- settings_DE$res_settings$nom_DE
            settings_DE$denom_DE <- settings_DE$res_settings$denom_DE
            settings_DE$batchVar1 <- settings_DE$res_settings$batchVar1
            settings_DE$batchVar2 <- settings_DE$res_settings$batchVar2
            settings_DE$IRmode_DE <- settings_DE$res_settings$IRmode_DE

            if("filters" %in% names(load_DE)) {
                settings_DE$filters <- load_DE$filters
            }
        })

        observeEvent(input$clear_selected_DE, {
            req(settings_DE$res)
            req(input$DT_DE_rows_selected)
            settings_DE$command_selected <- NULL
            DT::dataTableProxy("DT_DE") %>% DT::selectRows(NULL)
            settings_DE$DT_DE_rows_selected <- NULL
        })
        
        observeEvent(get_selected_diag(), {
            settings_DE$command_selected <- get_selected_diag()
        })
        
        observeEvent(get_selected_volc(), {
            settings_DE$command_selected <- get_selected_volc()
        })
        
        # Double check command_selected is refreshed in the table
        observeEvent(settings_DE$command_selected, {
            if(
                    is_valid(settings_DE$command_selected) && 
                    !identical(settings_DE$command_selected, 
                        input$DT_DE_rows_selected)
            ) {
                DT::dataTableProxy("DT_DE") %>% 
                    DT::selectRows(settings_DE$command_selected)        
            } else if(
                !identical(settings_DE$command_selected,
                    input$DT_DE_rows_selected)
            ) {
                DT::dataTableProxy("DT_DE") %>% DT::selectRows(NULL)
            }
        })
    
        return(settings_DE)
    })
}

.get_discrete_cats <- function(se) {
    if(!is(se, "NxtSE")) return(NULL)
    colData <- colData(se)
    if(ncol(colData) == 0) return(NULL)
    
    ret <- c()
    for(colcat in colnames(colData)) {
        if(is(colData[, colcat], "factor")) {
            ret <- c(ret, colcat)
        } else if(is(colData[, colcat], "numeric")) {
            ret <- c(ret, colcat)
        }
    }
    return(ret)
}