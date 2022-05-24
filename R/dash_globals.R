.ui_notice <- function(top = 50, left = 50, maxwidth = 250) {
    tags$head(tags$style(
        paste0(
            "#shiny-notification-panel {",
            "position: fixed;",
            "top: ", as.character(top), "%;", "bottom: unset;",
            "left: ", as.character(left), "%;", "right: unset;",
            "align-items: center;",
            "margin-left: auto;",
            "margin-right: auto;",
            "width: 100%;",
            "max-width: ", as.character(maxwidth), "px;}"
        )
    ))
}
