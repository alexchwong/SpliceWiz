covPlotly <- function(
        fig = list(),
        covTrack = list(),
        diffTrack = list(),
        annoTrack = list(),
        vLayout = c(6,1,2)
) {
    obj <- new("covPlotly",
        fig = fig,
        covTrack = covTrack,
        diffTrack = diffTrack,
        annoTrack = annoTrack,
        vLayout = vLayout
    )
    obj
}


setMethod("show", "covPlotly", function(object) {
    if(length(object@fig) == 0) return(NULL)
    if(!is(object@fig[[1]], "plotly")) return(NULL)
    p <- object@fig[[1]]
    
    for(i in seq_len(length(p$x$data))) {
        p$x$data[[i]]$hoveron <- NULL
    }
    
    show(p)
})
