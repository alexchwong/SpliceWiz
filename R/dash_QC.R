ui_qc <- function(id) {
    ns <- NS(id)
    wellPanel(
        .ui_notice(),
        fluidRow(
            column(2,
                selectInput(ns('QCmode'), 'Mode', width = '100%',
                    choices = c("Summary Table", "PCA", "Graphs")), 
                conditionalPanel(
                    ns = ns,
                    condition = "['Graphs'].indexOf(input.QCmode) >= 0",
                    selectInput(ns('QC_xaxis'), 'X axis', width = '100%',
                        choices = c("(none)")), 
                    selectInput(ns('QC_yaxis'), 'Y axis', width = '100%',
                        choices = c("(none)"))
                )
            ),
            column(10,
                conditionalPanel(
                    ns = ns,
                    condition = "['Summary Table'].indexOf(input.QCmode) >= 0",
                    div(style = 'overflow-x: scroll',  
                        DT::dataTableOutput(ns('DT_QC'))
                    )
                ),
                conditionalPanel(
                    ns = ns,
                    condition = "['Summary Table'].indexOf(input.QCmode) < 0",
                    plotlyOutput(ns("QC_plot"), height = "800px")
                )
            )
        )
    )
}

server_qc <- function(id, refresh_tab, get_se_path, get_df) {
    moduleServer(id, function(input, output, session) {
        settings_QC <- setreactive_QC()
        
        observeEvent(refresh_tab(), {
            updateSelectInput(session = session, inputId = "QC_xaxis",
                choices = "(none)")
            updateSelectInput(session = session, inputId = "QC_yaxis",
                choices = "(none)")
            output$DT_QC <- DT::renderDataTable(NULL)
            
            if(file.exists(file.path(get_se_path(), "stats.fst"))) {
                settings_QC$QC <- as.data.table(read.fst(
                    file.path(get_se_path(), "stats.fst")))               
                settings_QC$QC <- 
                    merge(as.data.table(get_df()), settings_QC$QC, all = TRUE)
            }
            output$DT_QC <- DT::renderDataTable({
                DT::datatable(
                    as.data.frame(settings_QC$QC),
                    class = 'cell-border stripe',
                    rownames = settings_QC$QC$sample,
                    filter = 'top'
                )
            })
            if(is_valid(settings_QC$QC)) {
                choices <- colnames(settings_QC$QC)
                choices <- choices[
                    !(choices %in% colnames(get_df()))
                ]
                choices <- choices[!(choices %in% 
                    c("paired", "strand", "path")
                )]
                choices <-c("(none)", choices)
                updateSelectInput(session = session, inputId = "QC_xaxis",
                    choices = choices)
                updateSelectInput(session = session, inputId = "QC_yaxis",
                    choices = choices)
            }      
        })
        
        observeEvent({ list(
            input$QCmode,
            input$QC_xaxis,
            input$QC_yaxis
        )}, {
            req(settings_QC$QC)
            choices <- colnames(settings_QC$QC)
            choices <- choices[!(choices %in% colnames(get_df()))]
            choices <- choices[!(choices %in% 
                c("paired", "strand", "path")
            )]
            output <- QC_update_plot(settings_QC$QC, choices, input$QCmode,
                input$QC_xaxis, input$QC_yaxis, output)
        })

    })
}

QC_update_plot <- function(QC, QC_cols, mode, x_axis, y_axis, output) {
    df <- as.data.frame(QC)
    rownames(df) <- df$sample
    if(mode == "PCA") {
        mat <- as.matrix(df[, QC_cols])
        rownames(mat) <- df$sample
        output$QC_plot <- renderPlotly({
            QC_PCA(mat)
        })
    } else if(mode == "Graphs") {
        output$QC_plot <- renderPlotly({
            validate(need(is_valid(x_axis) | is_valid(y_axis),
                "Specify X or Y axis"))
            if(is_valid(x_axis) & is_valid(y_axis)) {
                QC_Scatter_XY(df, x_axis, y_axis)
            } else if(is_valid(x_axis)) {
                QC_Bar_X(df, x_axis)
            } else if(is_valid(y_axis)) {
                QC_Bar_Y(df, y_axis)
            }
        })
    }    
    return(output)
}

QC_PCA <- function(mat) {
    colVar <- colVars(mat)
    mat <- mat[,colVar > 0]
    PCA <- prcomp(mat, scale. = TRUE)
    print(
        ggplotly(
            ggplot(as.data.frame(PCA$x), 
                aes(x = get("PC1"), y = get("PC2"), 
                    text = rownames(PCA$x))
            ) + geom_point() + 
            geom_text(aes(label = rownames(PCA$x))),
            tooltip = "text"
        )
    )
}

QC_Scatter_XY <- function(QC, x_axis, y_axis) {
    df.plot <- data.frame(
        sample = QC$sample,
        xaxis = unname(unlist(QC[,x_axis])),
        yaxis = unname(unlist(QC[,y_axis]))
    )
    colnames(df.plot)[2:3] <- c(
        x_axis, y_axis
    )
    print(ggplotly(
        ggplot(df.plot, aes_string(
            x = x_axis, y = y_axis,
            text = "sample")) +
        geom_point() + geom_text(aes(label = sample)),
        tooltip = "text"
    ))
}

QC_Bar_X <- function(QC, x_axis) {
    df.plot <- data.frame(sample = QC$sample,
        xaxis = unname(unlist(QC[,x_axis]))
    )
    colnames(df.plot)[2] <- x_axis
    print(ggplotly(
        ggplot(df.plot, aes_string(
            x = x_axis, y = "sample",
            text = "sample")) +
        geom_bar(stat="identity"),
        tooltip = "text"
    ))  
}

QC_Bar_Y <- function(QC, y_axis) {
    df.plot <- data.frame(sample = QC$sample,
        yaxis = unname(unlist(QC[,y_axis]))
    )
    colnames(df.plot)[2] <- y_axis
    print(ggplotly(
        ggplot(df.plot, aes_string(
            y = y_axis, x = "sample",
            text = "sample")) +
        geom_bar(stat="identity"),
        tooltip = "text"
    ))
}