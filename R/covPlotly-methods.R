#' Container for plotly-based coverage plots
#'
#' @param object A covPlotly object
#' @param resolution The number of horizontal "pixels" or data-points to plot.
#'   This is calculated per sub-plot. Smaller numbers lead to lower resolution
#'   but faster plots.
#' @return
#'   For `show()`: A plotly object synthesised by `plotView()`
#'   For `getExonRanges()`: A named GRanges object containing exon ranges
#'   For `showExons()`: A named GRanges object containing exon ranges, and
#'     additionally "shows" the plotly coverage plot with annotation replaced
#'     by named exons
#'   For `setResolution()` Returns the `covPlotly` object with addition of
#'     resolution set by the corresponding parameter. When `show()` is called,
#'     the plotly object with the new coverage resolution will be displayed.
#' @examples
#' se <- SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#'
#' # Assign annotation of the experimental conditions
#' colData(se)$treatment <- rep(c("A", "B"), each = 3)
#'
#' # Retrieve coverage data for all samples for the gene "SRSF3" (and surrounds)
#' 
#' dataObj <- getCoverageData(
#'     se,
#'     Gene = "SRSF3",
#'     tracks = colnames(se)
#' )
#'
#' plotObj_samples <- getPlotObject(
#'     dataObj,
#'     Event = "SE:SRSF3-203-exon4;SRSF3-202-int3"
#' )
#' 
#' if(interactive()) {
#'
#'     # Create covPlotly object by setting `usePlotly = TRUE`
#'     p <- plotView(plotObj_samples, usePlotly = TRUE)
#'
#'     # Display plotly plot
#'     show(p)
#'
#'     # Set resolution to 2000; display new plot
#'     p <- setResolution(p, resolution = 2000)
#'     show(p)
#'
#'     # Display exon annotation along with generated plot;
#'     # - also returns GRanges object
#'     gr <- showExons(p)
#' }
#' 
#' @name covPlotly-class
#' @aliases
#' getExonRanges getExonRanges,covPlotly-method
#' setResolution setResolution,covPlotly-method
#' showExons showExons,covPlotly-method
#' @seealso [plotView]
NULL

covPlotly <- function(
        fig = list(),
        args = list(),
        covTrack = list(),
        diffTrack = list(),
        annoTrack = list(),
        exonTrack = list(),
        vLayout = c(6,1,2)
) {
    obj <- new("covPlotly",
        fig = fig,
        args = args,
        covTrack = covTrack,
        diffTrack = diffTrack,
        annoTrack = annoTrack,
        exonTrack = exonTrack,
        vLayout = vLayout
    )
    obj
}

setMethod("show", "covPlotly", function(object) {
    if(length(object@fig) < 1) return(NULL)
    if(!is(object@fig[[1]], "plotly")) return(NULL)
    if(length(object@fig) == 1) {
        show(object@fig[[1]])
    } else {
        p <- object@fig[[2]]
        
        if("resolution" %in% names(object@args)) {
            resolution <- object@args[["resolution"]]
        } else {
            resolution <- 5000
        }
        
        if("xrange" %in% names(object@args)) {
            rangeStart <- min(object@args[["xrange"]])
            rangeEnd <- max(object@args[["xrange"]])
            rangeWidth <- rangeEnd - rangeStart

            okCoords <- .pV_getAllowedCoords(
                rangeStart - rangeWidth, rangeEnd + rangeWidth,
                object@args[["reservedCoords"]], resolution
            )
            
            # cull x-coordinates that are not reserved
            if("covTrackPos" %in% names(object@args)) {
                for(j in seq_len(length(object@args[["covTrackPos"]]))) {
                    curTrack <- object@args[["covTrackPos"]][[j]] - 1
                    nTracks <- object@args[["numCovTraces"]][[j]]
                    
                    for(k in seq(3, 2 + nTracks * 2)) {
                        DT <- data.table(
                            x = p$x$data[[curTrack + k]]$x,
                            y = p$x$data[[curTrack + k]]$y,
                            text = p$x$data[[curTrack + k]]$text
                        )
                        if(!is.null(okCoords)) DT <- DT[get("x") %in% okCoords]
                        p$x$data[[curTrack + k]]$x <- DT$x
                        p$x$data[[curTrack + k]]$y <- DT$y
                        p$x$data[[curTrack + k]]$text <- DT$text
                    }
                }
            }
            if("diffTrackPos" %in% names(object@args)) {
                for(j in seq_len(length(object@args[["diffTrackPos"]]))) {
                    curTrack <- object@args[["diffTrackPos"]][[j]] - 1
                    nTracks <- object@args[["numDiffTraces"]][[j]]
                    DT <- data.table(
                        x = p$x$data[[curTrack + 1]]$x,
                        y = p$x$data[[curTrack + 1]]$y,
                        text = p$x$data[[curTrack + 1]]$text
                    )
                    if(!is.null(okCoords)) DT <- DT[get("x") %in% okCoords]
                    p$x$data[[curTrack + 1]]$x <- DT$x
                    p$x$data[[curTrack + 1]]$y <- DT$y
                    p$x$data[[curTrack + 1]]$text <- DT$text
                }
            }
        }

        show(p)
    }
})

#' @describeIn covPlotly-class Returns a named GRanges object containing exon
#' ranges, without showing the associated plotly object
#' @export
setMethod("getExonRanges", "covPlotly", function(object) {
    return(object@args[["exonRanges"]])
})

#' @describeIn covPlotly-class Returns a covPlotly object after setting
#' the output resolution of the plotly-based coverage plots.
#' @param resolution How many horizontal pixels of resolution should be shown
#'   in the final plotly object. Set to `0` to disable.
#' @export
setMethod("setResolution", "covPlotly", function(object, resolution) {
    if(is.numeric(resolution)) {
        object@args[["resolution"]] <- resolution
    }
    return(object)
})

#' @describeIn covPlotly-class Returns a named GRanges object containing exon
#' ranges, and shows the plotly object with the annotation track showing the
#' named exons
#' @export
setMethod("showExons", "covPlotly", function(object) {
    if(length(object@fig) < 1) return(NULL)
    if(!is(object@fig[[1]], "plotly")) return(NULL)
    if(length(object@annoTrack) < 1) return(NULL)
    if(length(object@fig) == 1) {
        show(object@fig[[1]])
    } else {
        # Inject exon ranges
        p <- object

        p <- injectPlotData(p, p@args[["annoTrackPos"]], 
            p@exonTrack[[1]][["dataList"]],
            p@exonTrack[[1]][["layoutList"]][["xtitle"]]
        )

        fig <- p@fig[[2]]
        
        if("resolution" %in% names(p@args)) {
            resolution <- p@args[["resolution"]]
        } else {
            resolution <- 5000
        }
        
        if("xrange" %in% names(p@args)) {
            rangeStart <- min(p@args[["xrange"]])
            rangeEnd <- max(p@args[["xrange"]])
            rangeWidth <- rangeEnd - rangeStart
            okCoords <- round(seq(
                from = rangeStart - rangeWidth,
                to = rangeEnd + rangeWidth,
                length.out = resolution - length(p@args[["reservedCoords"]])
            ))
            okCoords <- sort(unique(c(okCoords, p@args[["reservedCoords"]])))
            
            # cull x-coordinates that are not reserved
            if("covTrackPos" %in% names(p@args)) {
                for(j in seq_len(length(p@args[["covTrackPos"]]))) {
                    curTrack <- p@args[["covTrackPos"]][[j]] - 1
                    nTracks <- p@args[["numCovTraces"]][[j]]
                    
                    for(k in seq(3, 2 + nTracks * 2)) {
                        DT <- data.table(
                            x = fig$x$data[[curTrack + k]]$x,
                            y = fig$x$data[[curTrack + k]]$y,
                            text = fig$x$data[[curTrack + k]]$text
                        )
                        DT <- DT[get("x") %in% okCoords]
                        fig$x$data[[curTrack + k]]$x <- DT$x
                        fig$x$data[[curTrack + k]]$y <- DT$y
                        fig$x$data[[curTrack + k]]$text <- DT$text
                    }
                }
            }
            if("diffTrackPos" %in% names(p@args)) {
                for(j in seq_len(length(p@args[["diffTrackPos"]]))) {
                    curTrack <- p@args[["diffTrackPos"]][[j]] - 1
                    nTracks <- p@args[["numDiffTraces"]][[j]]
                    DT <- data.table(
                        x = fig$x$data[[curTrack + 1]]$x,
                        y = fig$x$data[[curTrack + 1]]$y,
                        text = fig$x$data[[curTrack + 1]]$text
                    )
                    DT <- DT[get("x") %in% okCoords]
                    fig$x$data[[curTrack + 1]]$x <- DT$x
                    fig$x$data[[curTrack + 1]]$y <- DT$y
                    fig$x$data[[curTrack + 1]]$text <- DT$text
                }
            }
  
            show(fig)
            return(object@args[["exonRanges"]])
        }       
    }
})

# for coverage, and also for introns
addLineTrace <- function(fig, colorCode = "#000000") {
    if(!is(fig, "plotly")) return(fig)
    
    fig %>% add_trace(
        type = "scatter", mode = "lines",
        x = c(1,2), y = c(1,2), text = rep("test", 2),
        hoveron = "points", hoverinfo = I("text"),
        line = list(color = colorCode),
        visible = FALSE,
        showlegend = FALSE
    )
}

addRibbonTrace <- function(fig, colorCode = "#000000") {
    if(!is(fig, "plotly")) return(fig)
    
    fig %>% add_ribbons(
        type = "scatter", mode = "lines",
        x = c(1,2), ymin = c(1,2) - 0.2, ymax = c(1,2) + 0.2, 
        text = rep("test", 2),
        hoveron = "points", hoverinfo = I("text"),
        line = list(color = colorCode),
        color = I(colorCode), opacity = 0.2,
        visible = FALSE,
        showlegend = FALSE
    )
}

addJuncTrace <- function(fig, colorCode = "rgb(255, 100, 100)") {
    if(!is(fig, "plotly")) return(fig)
    
    fig %>% add_trace(
        type = "scatter", mode = "lines",
        x = c(1,2), y = c(1,2), text = rep("test", 2),
        hoveron = "points", hoverinfo = I("text"),
        line = list(color = colorCode, width = 0.5),
        visible = FALSE,
        showlegend = FALSE
    )
}

addTextTrace <- function(fig) {
    if(!is(fig, "plotly")) return(fig)
    
    fig %>% add_trace(
        type = "scatter", mode = "text", textposition = "middle",
        x = c(1,2), y = c(1,2), text = rep("test", 2),
        hoverinfo = I("text"),
        visible = FALSE,
        showlegend = FALSE
    )
}

addExonTrace <- function(fig, colorCode = "rgb(255, 100, 100)") {
    if(!is(fig, "plotly")) return(fig)
    
    fig %>% add_trace(
        x = c(1,1,2,2,1),
        y = c(1,2,2,1,1),
        text = rep("test", 5),
        type = 'scatter', mode = 'lines',
        hoveron = "points", hoverinfo = 'text',
        line = list(
            color = "transparent"
        ),
        fill = "toself",
        fillcolor = colorCode,
        visible = FALSE,
        showlegend = FALSE
    )
}

addCovTrack <- function(n_traces) {
    fig <- plot_ly() 
    cols <- scales::hue_pal()(n_traces)
    if(n_traces == 1) cols <- "#000000"
    
    # always junctions first
    fig <- fig %>%
        addJuncTrace() %>% addTextTrace()

    for(colorCode in cols) {
        fig <- fig %>%
            addLineTrace(colorCode) %>%
            addRibbonTrace(colorCode)
    }
    return(fig)
}


addDiffTrack <- function(n_traces) {
    fig <- plot_ly()
    cols <- scales::hue_pal()(n_traces)
    if(n_traces == 1) cols <- "#000000"
    
    for(colorCode in cols) {
        fig <- fig %>%
            addLineTrace(colorCode)
    }
    return(fig)
}

# traces
# line: black, blue, red, purple
# exons: black, blue, red, purple
addAnnoTrack <- function() {
    fig <- plot_ly()
    colors <- c(
        "rgba(0,0,0,1)", "rgba(0,0,255,1)",
        "rgba(255,0,0,1)", "rgba(255,0,255,1)"
    )
    for(col in colors) {
        fig <- fig %>% addLineTrace(col)
    }
    for(col in colors) {
        fig <- fig %>% addExonTrace(col)
    }
    fig <- fig %>% addTextTrace()
    return(fig)
}

knitPlotly <- function(
    p,
    numCovTraces = c(1,1),
    numDiffTraces = c(1),
    vLayout = c(6,1,2)
) {
    if(!is(p, "covPlotly")) return(p)

    if(
        identical(p@args[["numCovTraces"]], numCovTraces) &
        identical(p@args[["numDiffTraces"]], numDiffTraces) &
        identical(p@args[["vLayout"]], vLayout)
    ) {
        p@fig[[2]] <- p@fig[[1]]
        return(p)
    }

    figList <- list()
    figCount <- 0

    p@args[["vLayout"]] <- vLayout
    p@args[["numCovTraces"]] <- numCovTraces
    p@args[["numDiffTraces"]] <- numDiffTraces
    
    p@args[["covTrackPos"]] <- c()
    p@args[["diffTrackPos"]] <- c()
    p@args[["annoTrackPos"]] <- c()
    
    traceCount <- 0
    for(n in numCovTraces) {
        figCount <- figCount + 1
        figList[[figCount]] <- addCovTrack(n)
        p@args[["covTrackPos"]] <- c(p@args[["covTrackPos"]], traceCount + 1)
        traceCount <- traceCount + (2 + 2*n)
    }
    for(n in numDiffTraces) {
        figCount <- figCount + 1
        figList[[figCount]] <- addDiffTrack(n)
        p@args[["diffTrackPos"]] <- c(p@args[["diffTrackPos"]], traceCount + 1)
        traceCount <- traceCount + (n)
    }
    # Always one annotation track
    figCount <- figCount + 1
    figList[[figCount]] <- addAnnoTrack()
    p@args[["annoTrackPos"]] <- traceCount + 1

    vHeights <- vLayout
    vlNorm <- c()
    if(length(numCovTraces) > 0) {
        vlNorm <- c(vlNorm, 
            rep(
                vHeights[1] / length(numCovTraces), 
                length(numCovTraces)
            )
        )
    }
    if(length(numDiffTraces) > 0) {
        vlNorm <- c(vlNorm, 
            rep(
                vHeights[2] / length(numDiffTraces), 
                length(numDiffTraces)
            )
        )
    }
    vlNorm <- c(vlNorm, vHeights[3])
    vlNorm <- vlNorm / sum(vlNorm)
    
    fig <- subplot(
        figList, nrows = length(vlNorm),
        shareX = TRUE, titleY = TRUE,
        heights = vlNorm
    )
    p@fig[[1]] <- fig
    p@fig[[2]] <- p@fig[[1]]
    return(p)
}

injectPlotData <- function(p, trackPos, dataList, trackName = "trackN") {
    if(!is(p, "covPlotly")) return(p)
    if(length(p@fig) != 2) return(p)
    if(!is(p@fig[[2]], "plotly")) return(p)
    
    fig <- p@fig[[2]]   
    if(trackPos + length(dataList) - 1 > length(fig$x$data)) {
        .log("Data is longer than end of subplot")
    }

    curTrack <- trackPos - 1
    for(i in seq_len(length(dataList))) {
        curTrack <- curTrack + 1
        if(length(dataList[[i]]) > 0) {
            fig$x$data[[curTrack]]$x <- dataList[[i]]$x
            fig$x$data[[curTrack]]$y <- dataList[[i]]$y
            fig$x$data[[curTrack]]$text <- dataList[[i]]$text
            if("hovertemplate" %in% names(dataList[[i]])) {
                fig$x$data[[curTrack]]$hovertemplate <- 
                    dataList[[i]]$hovertemplate
                fig$x$data[[curTrack]]$hoverinfo <- NULL
            }
            fig$x$data[[curTrack]]$visible <- TRUE
            if("showlegend" %in% names(dataList[[i]])) {
                fig$x$data[[curTrack]]$showlegend <- 
                    dataList[[i]]$showlegend
                fig$x$data[[curTrack]]$name <- dataList[[i]]$name
            } else {
                fig$x$data[[curTrack]]$showlegend <- FALSE
                fig$x$data[[curTrack]]$name <- trackName
            }
            
            fig$x$data[[curTrack]]$hoverinfo <- "text"
        } else {
            fig$x$data[[curTrack]]$x <- c(1,2)
            fig$x$data[[curTrack]]$y <- c(1,2)
            fig$x$data[[curTrack]]$text <- rep("test", 2)
            fig$x$data[[curTrack]]$visible <- FALSE
            fig$x$data[[curTrack]]$showlegend <- FALSE
        }
    }
    p@fig[[2]] <- fig
    return(p)
}

adjustXrange <- function(p, rangeX) {
    if(!is(p, "covPlotly")) return(p)
    if(length(p@fig) != 2) return(p)
    if(!is(p@fig[[2]], "plotly")) return(p)
    
    p@fig[[2]]$x$layout$xaxis[["range"]] <- rangeX
    p@args[["xrange"]] <- rangeX
    return(p)
}

adjustXtitle <- function(p, titleName) {
    if(!is(p, "covPlotly")) return(p)
    if(length(p@fig) != 2) return(p)
    if(!is(p@fig[[2]], "plotly")) return(p)
    
    p@fig[[2]]$x$layout$xaxis[["title"]] <- titleName
    return(p)
}

adjustYrange <- function(p, rangeY, trackNum) {
    if(!is(p, "covPlotly")) return(p)
    if(length(p@fig) != 2) return(p)
    if(!is(p@fig[[2]], "plotly")) return(p)
    
    if(trackNum == 1) {
        axisName <- "yaxis"
    } else {
        axisName <- paste0("yaxis", as.character(trackNum))
    }
    
    p@fig[[2]]$x$layout[[axisName]][["range"]] <- rangeY
    return(p)
}

adjustYtitle <- function(p, titleName, trackNum) {
    if(!is(p, "covPlotly")) return(p)
    if(length(p@fig) != 2) return(p)
    if(!is(p@fig[[2]], "plotly")) return(p)

    if(trackNum == 1) {
        axisName <- "yaxis"
    } else {
        axisName <- paste0("yaxis", as.character(trackNum))
    }

    p@fig[[2]]$x$layout[[axisName]][["title"]] <- titleName
    return(p)
}

fixYrange <- function(p, trackNum) {
    if(!is(p, "covPlotly")) return(p)
    if(length(p@fig) != 2) return(p)
    if(!is(p@fig[[2]], "plotly")) return(p)

    if(trackNum == 1) {
        axisName <- "yaxis"
    } else {
        axisName <- paste0("yaxis", as.character(trackNum))
    }

    p@fig[[2]]$x$layout[[axisName]][["fixedrange"]] <- TRUE
    return(p)
}
