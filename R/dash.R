#' Launches the SpliceWiz Graphics User Interface using Shiny Dashboard
#' 
#' This function launches the SpliceWiz interactive app using Shiny Dashboard
#' @param mode (default `"dialog"`) `"dialog"` displays SpliceWiz in a dialog 
#'   box with specified width and height. `"browser"` opens SpliceWiz in a 
#'   browser-like resizable window.
#' @param res (default "1080p") Sets width and height of the app to pre-defined
#'   dimensions. Possible options are "720p, "960p", "1080p", "1440p", which
#'   specifies the height of the app. All are displayed in aspect ratio 16x9
#' @return Runs an interactive shinydashboard SpliceWiz app with the specified
#'   mode.
#' @examples
#' if(interactive()) {
#' 
#' # Launches interactive ShinyDashboard SpliceWiz app as fixed-size dialog box
#' # 1080p = 1920 x 1080 pixels
#'     spliceWiz(mode = "dialog", res = "1080p") 
#'
#' # Launches interactive ShinyDashboard SpliceWiz app as browser window
#'     spliceWiz(mode = "browser")
#' 
#' }
#' 
#' @md
#' @export
spliceWiz <- function(
        mode = c("dialog", "browser"),
        res = c("1080p", "720p", "960p", "1440p")
) {
    if(!interactive()) {
        .log(paste("In spliceWiz(), SpliceWiz App",
            "can only be run in interactive mode (i.e. RStudio)."))
    }
    mode = match.arg(mode)
    res = match.arg(res)
    
    ui_dash <- dashboardPage(
        dashboardHeader(title = "SpliceWiz"),
        ui_sidebar(),
        dashboardBody(
            tabItems(
                ui_tab_title(),     # "navTitle"

                ui_tab_ref_new(),   # "navRef_New"    

                ui_tab_expr(),      # "navExpr"
                
                ui_tab_expr_load(), # "navExprLoad"
                ui_tab_qc(),        # "navQC"
                ui_tab_filter(),    # "navFilter"
                ui_tab_analyse(),   # "navAnalyse"

                ui_tab_volcano(),   # "navVolcano"
                ui_tab_diag(),      # "navDiag"
                ui_tab_heatmap(),   # "navHeatmap"
                ui_tab_coverage()   # "navCoverage"
            )
        )
    )
    height = 1080
    width = 1920
    if(res == "720p") {
        height = 720
        width = 1280
    } else if(res == "960p") {
        height = 960
        width = 1700    
    } else if(res == "1440p") {
        height = 1440
        width = 2560
    }
    
    if(mode == "browser") {
        runApp(shinyApp(ui_dash, dash_server))
    } else {
        runGadget(
            shinyApp(ui_dash, dash_server),
            viewer = dialogViewer('SpliceWiz', 
                width = width, height = height)
        )
    }
}