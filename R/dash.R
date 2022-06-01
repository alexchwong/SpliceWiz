#' Launches the SpliceWiz Graphics User Interface (GUI) using Shiny Dashboard
#' 
#' This function launches the SpliceWiz interactive app using Shiny Dashboard
#' This is (by default) a dialog window within the RStudio application with
#' the resolution specified by the `res` parameter. Alternatively, setting
#' `mode = "browser"` will launch a resizable browser window (using the default
#' internet browser). The demo mode can be launched by setting `demo = TRUE`.
#' See the [SpliceWiz Quick-Start](../doc/SW_QuickStart.html) for a guide
#' to using the SpliceWiz GUI. 
#' @param mode (default `"dialog"`) `"dialog"` displays SpliceWiz in a dialog 
#'   box with specified width and height. `"browser"` opens SpliceWiz in a 
#'   browser-like resizable window.
#' @param res (default "1080p") Sets width and height of the app to pre-defined
#'   dimensions. Possible options are "720p, "960p", "1080p", "1440p", which
#'   specifies the height of the app. All are displayed in aspect ratio 16x9
#' @param demo (default FALSE) If set to `TRUE`, SpliceWiz will place demo
#'   reference and BAM files into the temporary directory.
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
#' @name Graphics-User-Interface
#' @aliases
#' GUI spliceWiz
#' @md
NULL

#' @describeIn Graphics-User-Interface
#' @export
spliceWiz <- function(
        mode = c("dialog", "browser"),
        res = c("1080p", "720p", "960p", "1440p"),
        demo = FALSE
) {
    if(!interactive()) {
        .log(paste("In spliceWiz(), SpliceWiz App",
            "can only be run in interactive mode (i.e. RStudio)."))
    }
    mode <- match.arg(mode)
    res <- match.arg(res)
    
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
    height <- 1080
    width <- 1920
    if(res == "720p") {
        height <- 720
        width <- 1280
    } else if(res == "960p") {
        height <- 960
        width <- 1700    
    } else if(res == "1440p") {
        height <- 1440
        width <- 2560
    }
    
    if(demo) {
        if(!dir.exists(file.path(tempdir(), "Reference")))
            dir.create(file.path(tempdir(), "Reference"))
        if(!dir.exists(file.path(tempdir(), "bams")))
            dir.create(file.path(tempdir(), "bams"))
        if(!dir.exists(file.path(tempdir(), "pb_output")))
            dir.create(file.path(tempdir(), "pb_output"))
        if(!dir.exists(file.path(tempdir(), "NxtSE")))
            dir.create(file.path(tempdir(), "NxtSE"))
        ret <- example_bams(path = file.path(tempdir(), "bams"))
        
        if(is.null(ret)) {
            .log("Error creating demo BAM files", "message")
        } else {
            .log("Reference and BAM files placed into temporary directory", 
                "message")
            fwrite(
                data.frame(
                    sample = c("02H003", "02H025", "02H026",
                        "02H033", "02H043", "02H046"),
                    condition = rep(c("A", "B"), each = 3),
                    batch = rep(c("K", "L", "M"), 2)
                ),
                file.path(tempdir(), "demo_annotations.csv")
            )        
            setwd(tempdir())
        }
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