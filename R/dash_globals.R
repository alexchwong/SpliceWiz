.ui_notice <- function(top = 50, left = 50, maxwidth = 250) {
    tags$head(tags$style(
        paste0(
            ".shiny-notification {",
            "position: fixed;",
            "top: ", as.character(top), "%;",
            "left: ", as.character(left), "%;",
            "align-items: center;",
            "margin-left: auto;",
            "margin-right: auto;",
            "width: 100%;",
            "max-width: ", as.character(maxwidth), "px;}"
        )
    ))
}