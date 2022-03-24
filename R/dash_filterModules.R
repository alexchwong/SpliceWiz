filterModule_UI <- function(id, label = "Counter") {
    ns <- NS(id)
    wellPanel(
        h5(label),  # e.g. "Filter #1"
        selectInput(ns("filterClass"), "Filter Class", 
            width = '100%', choices = c("(none)", "Annotation", "Data")),
        selectInput(ns("filterType"), "Filter Type", 
            width = '100%', choices = c("(none)")),
        conditionalPanel(ns = ns,
        condition = paste0("['TSL'].",
            "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_TSL_min"), 
                "TSL Threshold", 
                choices = seq_len(5), selected = 1)
        ),
        conditionalPanel(ns = ns,
            condition = paste0("['Consistency'].",
                "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_cons_max"), 
                "log-fold maximum", choices = seq(0.2, 5, by = 0.2), 
                selected = 1)
        ),
        conditionalPanel(ns = ns,
            condition = "['Coverage'].indexOf(input.filterType) >= 0",
            sliderInput(ns("slider_cov_min"), "Percent Coverage", 
                min = 0, max = 100, value = 80)
        ),
        conditionalPanel(ns = ns,
            condition = "['Depth'].indexOf(input.filterType) >= 0",
            shinyWidgets::sliderTextInput(ns("slider_depth_min"), 
                "Minimum", choices = c(1,2,3,5,10,20,30,50,100,200,300,500), 
                selected = 20),
        ),
        conditionalPanel(ns = ns,
            condition = "['Depth', 'Coverage'].indexOf(input.filterType) >= 0",
            tagList(
                shinyWidgets::sliderTextInput(ns("slider_mincond"), 
                    "Minimum Conditions Satisfy Criteria", 
                    choices = c(as.character(seq_len(8)), "All"), 
                    selected = "All"),
                selectInput(ns("select_conds"), "Condition", width = '100%',
                    choices = c("(none)")),
                sliderInput(ns("slider_pcTRUE"), 
                    "Percent samples per condition satisfying criteria", 
                    min = 0, max = 100, value = 80)
            )
        ),
        conditionalPanel(ns = ns,
            condition = paste0("['Coverage', 'Consistency'].",
                "indexOf(input.filterType) >= 0"),
            shinyWidgets::sliderTextInput(ns("slider_minDepth"), 
                "Signal Threshold to apply criteria", 
                choices = c(1,2,3,5,10,20,30,50,100,200,300,500), 
                selected = 20),
        ),
        conditionalPanel(ns = ns,
            condition = "['(none)'].indexOf(input.filterClass) < 0",
            selectInput(ns("EventType"), "Splice Type", width = '100%', 
                multiple = TRUE,
                choices = c("IR", "MXE", "SE", "A5SS", "A3SS",
                    "AFE", "ALE", "RI"))
        )
    )
}

filterModule_server <- function(id, filterdata, conditionList) {
    moduleServer(id, function(input, output, session) {
        #final <- reactiveValues(default = ASEFilter()) # initialize to defaults
        final <- reactiveVal(
            value = ASEFilter(
                filterClass = "(none)"
            )
        )
        # Observe whether colData of NxtSE changes
        
        observeEvent(conditionList(), {
            fCond <- final()@condition
            choices_conds <- c("(none)", conditionList())
            if(
                    # Valid condition
                    length(choices_conds) > 1 && is_valid(fCond) && 
                    fCond %in% choices_conds[-1]
            ) {
                updateSelectInput(
                    session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = fCond
                )
            } else if(
                is_valid(fCond) && 
                !(fCond %in% choices_conds)                
            ){
                # If condition is valid but not in column, reset it and return
                updateSelectInput(
                    session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)"
                )
                return(final)
            } else {
                updateSelectInput(
                    session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)"
                )            
            }
        })

        # inputs from final -> UI
        observeEvent(filterdata(), {
            if(is(filterdata(), "ASEFilter")) final(filterdata())

            class_choices <- c("(none)", "Annotation", "Data")
            type_choices <- c("(none)")

            fClass <- final()@filterClass
            if(is_valid(fClass) && fClass %in% class_choices) {
                if(fClass == "Annotation") {
                    type_choices <- c("Protein_Coding", "NMD", "TSL", 
                        "Terminus", "ExclusiveMXE")
                } else if(fClass == "Data") {
                    type_choices <- c("Depth", "Coverage", "Consistency")
                }
                updateSelectInput(session = session, 
                    inputId = "filterClass", choices = class_choices, 
                    selected = fClass)
            } else {
                # fClass == "" | fClass == "(none)"
                updateSelectInput(session = session, 
                    inputId = "filterClass", choices = class_choices)
                updateSelectInput(session = session, 
                    inputId = "filterType", choices = type_choices)
                # return(final)
            }
            
            fType <- final()@filterType
            if(is_valid(fType) && fType %in% type_choices) {
                updateSelectInput(session = session, inputId = "filterType", 
                    choices = type_choices, selected = fType)
            } else if(is_valid(fClass) && fClass %in% class_choices) {
                # fClass != "" & fClass != "(none)"
                updateSelectInput(session = session, inputId = "filterType", 
                    choices = type_choices) # Sets default fType if not set
                # return(final)
            } else {
                # Invalid fClass
                updateSelectInput(session = session, 
                    inputId = "filterClass", choices = class_choices)
                updateSelectInput(session = session, 
                    inputId = "filterType", choices = type_choices)
                # return(final)
            }

            feType <- final()@EventTypes
            eOptions <- c("IR", "MXE", "SE", "A3SS", "A5SS", "ALE", "AFE", "RI")
            
            # make sure feType is always valid
            if(length(feType) > 0) feType <- feType[feType %in% eOptions]
            if(length(feType) == 0) feType <- eOptions
            updateSelectInput(session = session, 
                inputId = "EventType", 
                selected = feType)

            fMin <- final()@minimum # always valid
            if(fType == "Depth") {
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_depth_min", 
                    selected = fMin)
            } else if(final()@filterType == "Coverage"){
                updateSliderInput(session = session, 
                    inputId = "slider_cov_min", 
                    value = fMin)
            } else  if(final()@filterType == "TSL"){
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_TSL_min", 
                    selected = fMin)
            }
            
            fMax <- final()@maximum # always valid
            shinyWidgets::updateSliderTextInput(
                session = session, inputId = "slider_cons_max", 
                selected = fMax)

            fmDepth <- final()@minDepth # always valid
            updateSelectInput(session = session, 
                inputId = "slider_minDepth", 
                selected = fmDepth)
            
            fmCond <- final()@minCond # always valid
            shinyWidgets::updateSliderTextInput(
                session = session, inputId = "slider_mincond", 
                selected = fmCond)
                
            choices_conds = c("(none)", conditionList())
            fCond <- final()@condition
            if(is_valid(fCond) && fCond %in% choices_conds) {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = fCond)
            } else {
                updateSelectInput(session = session, 
                    inputId = "select_conds", 
                    choices = choices_conds, 
                    selected = "(none)")            
            }
            
            fpcTRUE <- final()@pcTRUE
            updateSliderInput(session = session, 
                inputId = "slider_pcTRUE", 
                value = fpcTRUE)
        })

        # outputs from UI -> final
        observeEvent(input$filterClass, {
            obj <- final()
            obj@filterClass <- input$filterClass
            if(input$filterClass == "Annotation") {
                type_choices <- c("Protein_Coding", "NMD", "TSL", 
                    "Terminus", "ExclusiveMXE")
            } else if(input$filterClass == "Data") {
                type_choices <- c("Depth", "Coverage", "Consistency")
            } else {
                type_choices <- "(none)"
            }
            cur_choice <- obj@filterType
            if(is_valid(cur_choice) && cur_choice %in% type_choices) {
                updateSelectInput(session = session, 
                    inputId = "filterType", 
                    choices = type_choices, selected = cur_choice)
            } else {
                obj@filterType <- type_choices[1]
                updateSelectInput(session = session, 
                    inputId = "filterType", 
                    choices = type_choices)
            }
            final(obj)
        })
        observeEvent(input$filterType, {
            # final$trigger = NULL
            req(input$filterType)
            fType <- input$filterType
            
            obj <- final()
            obj@filterType <- fType

            fMin <- obj@minimum
            if(fType == "Depth") {
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_depth_min", 
                    selected = fMin)
            } else if(fType == "Coverage"){
                updateSliderInput(session = session, 
                    inputId = "slider_cov_min", 
                    value = fMin)
            } else if(fType == "TSL"){
                shinyWidgets::updateSliderTextInput(
                    session = session, inputId = "slider_TSL_min", 
                    selected = fMin)
            }
            final(obj)
        })
        observeEvent(input$slider_depth_min, {
            obj <- final()
            if(obj@filterType == "Depth") {
                obj@minimum = as.numeric(input$slider_depth_min)
            }
            final(obj)
        })
        observeEvent(input$slider_cov_min, {
            obj <- final()
            if(obj@filterType == "Coverage"){
                obj@minimum = as.numeric(input$slider_cov_min)
            }
            final(obj)
        })
        observeEvent(input$slider_TSL_min,{
            obj <- final()
            if(obj@filterType == "TSL"){
                obj@minimum = as.numeric(input$slider_TSL_min)
            }
            final(obj)
        })
        observeEvent(input$slider_cons_max,{
            obj <- final()
            obj@maximum = as.numeric(input$slider_cons_max)
            final(obj)
        })
        observeEvent(input$slider_minDepth,{
            obj <- final()
            obj@minDepth = as.numeric(input$slider_minDepth)
            final(obj)
        })
        observeEvent(input$slider_mincond,{
            obj <- final()
            if(input$slider_mincond == "All") {
                obj@minCond = -1
            } else {
                obj@minCond = as.numeric(input$slider_mincond)            
            }
            final(obj)
        })
        observeEvent(input$select_conds,{
            obj <- final()
            obj@condition = input$select_conds
            final(obj)
        })
        observeEvent(input$slider_pcTRUE,{
            obj <- final()
            obj@pcTRUE = as.numeric(input$slider_pcTRUE)
            final(obj)
        })
        observeEvent(input$EventType,{
            obj <- final()
            obj@EventTypes = input$EventType
            final(obj)
        })

        # Returns filter list from module
        return(final)
    })
}