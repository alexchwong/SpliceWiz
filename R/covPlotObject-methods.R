#' @describeIn covPlotObject-class Constructs a covPlotObject object.
#' For internal use only
#' @export
covPlotObject <- function(
        args = list(),
        cov = list(),
        junc = list(),
        cov_stats = list(),
        diff_stats = data.frame(),
        annotation = list()
) {
    obj <- new("covPlotObject",
        args = args,
        cov = cov,
        junc = junc,
        cov_stats = cov_stats,
        diff_stats = diff_stats,
        annotation = annotation
    )
    obj
}

#' @describeIn covPlotObject-class Creates a coverage plot using the stored
#'   data in the covPlotObject
#' @export
setMethod("plotCPO", c(x = "covPlotObject"), function(
    x, plotRanges = NULL,
    tracks = NULL,
    plotJunctions = TRUE,
    plotAnnotations = TRUE,
    plotAnnoSubTrack = TRUE,
    stackTracks = FALSE,
    plotDiffTrack = TRUE,
    interactive = FALSE,
    resolution = 10000,
    junctionThreshold = 0.01,
    ...
) {
    # Plot single track
    if(is.null(plotRanges)) {
        plotRanges <- GRanges(
            x@args[["view_chr"]],
            IRanges(
                x@args[["view_start"]], x@args[["view_end"]]
            )
        ) # strand is ignored
    }
    
    if(!is(plotRanges, "GRanges")) .log(paste(
        "In plot() for class covPlotObject,",
        "`plotRanges` must be a GRanges object, or NULL"
    ))

    if(
        !is.null(plotRanges) && is(plotRanges, "GRanges") &&
        any(as.character(seqnames(plotRanges)) != x@args[["view_chr"]])
    ) .log(paste(
        "In plot() for class covPlotObject,",
        "Some elements in `plotRanges` have seqnames that do not match that",
        "of covPlotObject"
    ))

    if(is.null(tracks)) tracks <- x@args[["tracks"]]

    if(!("tracks" %in% names(x@args))) .log(paste(
        "In plot() for class covPlotObject,",
        "attempted to specify `tracks` when there are none in covPlotObject"
    ))

    # Allot `tracks` to be index of tracks to use
    if(
            all(is.numeric(tracks)) && (
                all(tracks > 0) & all(tracks <= length(x@args[["tracks"]]))
            )
    ) {
        tracks <- unique(tracks)
        tracks <- x@args[["tracks"]][tracks]
    }

    if(!all(tracks %in% x@args[["tracks"]])) .log(paste(
        "In plot() for class covPlotObject,",
        "some `tracks` are not found in covPlotObject"
    ))
    
    if(!is.numeric(resolution) || resolution < 1000) {
        resolution <- 1000
    }
    
    # Structure of plot
    # | 1a || 1b || 1c |
    # | 2a || 2b || 2c |
    # | 3a || 3b || 3c |
    # | pa || pb || pc | - diff track
    # | A  || B  || C  | - annotation subtrack
    # | annotation     | - annotation full track

    covTrack <- list() # nested list
    diffTrack <- list()
    annoSubTrack <- list()
    annoFullTrack <- list() # list of 1
    
    subResolution <- resolution / ceiling(length(plotRanges))
    
    if(stackTracks) {
        for(j in seq_len(plotRanges)) {
            range_gr <- plotRanges[j]
            covTrack[[j]] <- .cPO_plotCoverage_multi(
                x, tracks, 
                range_gr = range_gr, 
                resolution = subResolution, 
                interactive = interactive
            )
        }
    } else {
        for(i in seq_len(length(tracks))) {
            covTrack[[i]] <- list()
            track <- tracks[i]
            for(j in seq_len(plotRanges)) {
                range_gr <- plotRanges[j]
                covTrack[[i]][[j]] <- .cPO_plotCoverage(
                    x, track, 
                    range_gr = range_gr, 
                    resolution = subResolution, 
                    jnThreshold = junctionThreshold,
                    plotJunctions = plotJunctions,
                    interactive = interactive
                )
            }
        }
    
    for(j in seq_len(plotRanges)) {
        if(stackTracks) {
        
        } else {
            for(i in seq_len(length(tracks))) {
                track <- tracks[i]
                covTrack[[i]][[j]] <- .cPO_plotCoverage(
                    x, track, 
                    range_gr = range_gr, 
                    resolution = subResolution, 
                    jnThreshold = junctionThreshold,
                    plotJunctions = plotJunctions,
                    interactive = interactive
                )
            }
        }
    }
})

.cPO_plotCoverage <- function(
    x, track, 
    range_gr, resolution, jnThreshold,
    plotJunctions = TRUE, interactive = FALSE
) {
    # plot mean coverage or raw coverage?
    plotMeanCov <- (length(x@cov_stats) > 0)

    # junction processing is identical
    juncAll <- x@junc
    junc <- x@junc[[track]]
    OL <- findOverlaps(junc, range_gr)
    if(length(from(OL)) > 0) {
        junc <- junc[unique(from(OL))]
    } else {
        junc <- NULL
    }

    if(plotMeanCov) {
        trackName <- paste(x@args[["condition"]], track)
        df <- x@cov_stats[[track]]
        df_sub <- df[,seq_len(2)] # only 2 columns
        if(nrow(df_sub) > 2 * resolution) {
            binwidth <- floor(nrow(df_sub) / resolution)
            df_sub <- .gCD_binCoverage(df_sub, binwidth, juncAll)
            
            # subset df by df_sub
            df <- df[df$x %in% round(df_sub$x),]
        }

        dfJn <- .cPO_jn_arcs(
            junc,
            arcHeight = 0.1 * max(df$mean),
            junctionThreshold = jnThreshold * max(df$mean)
        )
    
        if(x@args$ribbon_mode %in% c("ci", "sd", "sem")) {
            df$info <- paste(
                paste0("Coordinate: ", df$x), 
                paste0("Norm-Depth (mean): ", round(df$mean, 4)),
                paste0("Norm-Depth (", x@args$ribbon_mode, "): ", 
                    round(df$var, 4)),
                sep = "\n"
            )
        } else {
            df$info <- paste(
                paste0("Coordinate: ", df$x), 
                paste0("Norm-Depth (mean): ", round(df$mean, 4)),
                sep = "\n"
            )        
        }
    } else {
        trackName <- track
        
        df <- x@cov[[track]]
        if(nrow(df) > 2 * resolution) {
            binwidth <- floor(nrow(df) / resolution)
            df <- .gCD_binCoverage(df, binwidth, juncAll)
        }
        df$info <- paste(
            paste0("Coordinate: ", df$x), 
            paste0("Depth: ", df$depth),
            sep = "\n"
        )

        dfJn <- .cPO_jn_arcs(
            junc, 
            arcHeight = 0.1 * max(df$depth),
            junctionThreshold = jnThreshold * max(df$depth)
        )
    }
    df$group <- trackName

    # junc processing - common
    if(plotJunctions) {
        dtJn <- as.data.table(dfJn)
        dtJn <- dtJn[
            get("x") >= start(range_gr) &
            get("x") <= end(range_gr)
        ]
        dtJn[, c("xlabel", "ylabel") := list(
            mean(get("x")), mean(get("y"))), 
            by = "info"
        ]
        dtJn <- unique(dtJn[, 
            c("info", "value", "xlabel", "ylabel"), with = FALSE])
        dfJnSum <- as.data.frame(dtJn)                
    } else {
        dfJnSum <- NA
    }

    # plotting - common
        
    suppressWarnings({
        p <- ggplot(df, aes(text = get("info"))) +
            geom_hline(yintercept = 0)            
    })

    # plot ribbon (group cov only)
    if(plotMeanCov && x@args$ribbon_mode %in% c("ci", "sd", "sem")) {
        p <- p +
            geom_ribbon(data = df, alpha = 0.2,
                aes(
                    x = get("x"), y = get("mean"),
                    ymin = get("mean") - get("var"),
                    ymax = get("mean") + get("var"),
                    group = get("group")
                )
            )
    }
    
    # plot line
    if(plotMeanCov) {
        p <- p + geom_line(aes(
            x = get("x"), y = get("mean"), group = get("group")
        ))
    } else {
        p <- p + geom_line(aes(
            x = get("x"), y = get("depth"), group = get("group")
        ))
    }

    # plot y axis and format
    p <- p + theme_white_legend + labs(y = trackName)

    # plot junctions
    if(plotJunctions) {
        p <- p +
            geom_line(
                data = dfJn, 
                aes(x = get("x"), y = get("yarc"), group = get("info")), 
                color = "darkred"
            ) +
            geom_text(
                data = dfJnSum, 
                aes(
                    x = get("xlabel"), y = get("ylabel"), label = get("value")
                )
            )
        yrange <- c(0, 1.05 * max(
            c(layer_scales(p)$y$range$range[2], dfJn$yarc)
        ))
    } else {
        yrange <- c(0, 1.05 * layer_scales(p)$y$range$range[2])
    }
    
    pl <- NULL

        
    if(interactive) {
        pl <- ggplotly(p, tooltip = "text") %>% 
        layout(
            yaxis = list(
                range = yrange,
                rangemode = "tozero", 
                fixedrange = TRUE
            )
        )
        pl$x$data[[2]]$showlegend <- FALSE
        pl$x$data[[2]]$name <- trackName

        if(plotMeanCov) {
            pl$x$data[[3]]$showlegend <- FALSE
            pl$x$data[[3]]$name <- trackName        
        }

        return(pl)
    } else {
        p <- p +
            theme(axis.title.x = element_blank()) +
            labs(x = "", y = trackName)

        if(x@args[["reverseGenomeCoords"]]) {
            plotViewStart <- end(range_gr)
            plotViewEnd <- start(range_gr)
        } else {
            plotViewStart <- start(range_gr)
            plotViewEnd <- end(range_gr)            
        }

        p <- p + coord_cartesian(
            xlim = c(plotViewStart, plotViewEnd),
            ylim = yrange,
            expand = FALSE
        )     
        
        return(p)
    }
}

# Multiple traces on 1 track - junction plotting is not supported
.cPO_plotCoverage_multi <- function(
    x, tracks, trackName,
    range_gr, resolution, interactive = FALSE
) {
    # plot mean coverage or raw coverage?
    plotMeanCov <- (
        length(x@cov_stats) > 0 && all(tracks %in% names(x@cov_stats))
    )

    # junction processing is identical
    juncAll <- x@junc

    if(missing(trackName)) {
        trackName <- "Normalized Coverage"
    }

    df_all <- c()
    for(track in tracks) {
        if(plotMeanCov) {
            groupName <- paste(x@args[["condition"]], track)
            df <- x@cov_stats[[track]]
            df_sub <- df[,seq_len(2)] # only 2 columns
            if(nrow(df_sub) > 2 * resolution) {
                binwidth <- floor(nrow(df_sub) / resolution)
                df_sub <- .gCD_binCoverage(df_sub, binwidth, juncAll)
                
                # subset df by df_sub
                df <- df[df$x %in% round(df_sub$x),]
            }

            if(x@args$ribbon_mode %in% c("ci", "sd", "sem")) {
                df$info <- paste(
                    paste0(groupName, ":"),
                    paste0("Coordinate: ", df$x), 
                    paste0("Norm-Depth (mean): ", round(df$mean, 4)),
                    paste0("Norm-Depth (", x@args$ribbon_mode, "): ", 
                        round(df$var, 4)),
                    sep = "\n"
                )
            } else {
                df$info <- paste(
                    paste0("Coordinate: ", df$x), 
                    paste0("Norm-Depth (mean): ", round(df$mean, 4)),
                    sep = "\n"
                )        
            }
        } else {
            groupName <- track
            
            df <- x@cov[[track]]
            if(nrow(df) > 2 * resolution) {
                binwidth <- floor(nrow(df) / resolution)
                df <- .gCD_binCoverage(df, binwidth, juncAll)
            }
            df$info <- paste(
                paste0(groupName, ":"),
                paste0("Coordinate: ", df$x), 
                paste0("Depth: ", df$depth),
                sep = "\n"
            )
        }
        df$group <- groupName
        df_all <- rbind(df_all, df)
    }

    # plotting - common        
    suppressWarnings({
        p <- ggplot(df_all, aes(text = get("info"))) +
            geom_hline(yintercept = 0)            
    })

    # plot ribbon (group cov only)
    if(plotMeanCov && x@args$ribbon_mode %in% c("ci", "sd", "sem")) {
        p <- p +
            geom_ribbon(data = df_all, alpha = 0.2,
                aes(
                    x = get("x"), y = get("mean"),
                    ymin = get("mean") - get("var"),
                    ymax = get("mean") + get("var"),
                    group = get("group"),
                    fill = get("group")
                )
            )
    }
    
    # plot line
    if(plotMeanCov) {
        p <- p + geom_line(aes(
            x = get("x"), y = get("mean"), 
            group = get("group"), color = get("group")
        ))
    } else {
        p <- p + geom_line(aes(
            x = get("x"), y = get("depth"), 
            group = get("group"), color = get("group")
        ))
    }

    # plot y axis and format
    p <- p + theme_white_legend + labs(y = trackName)

    yrange <- c(0, 1.05 * layer_scales(p)$y$range$range[2])
    
    pl <- NULL
    if(interactive) {
        pl <- ggplotly(p, tooltip = "text") %>% 
        layout(
            yaxis = list(
                range = yrange,
                rangemode = "tozero", 
                fixedrange = TRUE
            )
        )
        pl$x$data[[2]]$showlegend <- FALSE
        pl$x$data[[2]]$name <- trackName

        if(plotMeanCov) {
            pl$x$data[[3]]$showlegend <- FALSE
            pl$x$data[[3]]$name <- trackName        
        }

        return(pl)
    } else {
        p <- p +
            theme(axis.title.x = element_blank()) +
            labs(x = "", y = trackName)

        if(x@args[["reverseGenomeCoords"]]) {
            plotViewStart <- end(range_gr)
            plotViewEnd <- start(range_gr)
        } else {
            plotViewStart <- start(range_gr)
            plotViewEnd <- end(range_gr)            
        }

        p <- p + coord_cartesian(
            xlim = c(plotViewStart, plotViewEnd),
            ylim = yrange,
            expand = FALSE
        )     
        
        return(p)
    }
}


.cPO_jn_arcs <- function(
        junc, 
        arcHeight = 0,
        junctionThreshold = 0
) {
    if(is.null(junc)) return(NA)
    
    df <- data.frame(
        juncName = paste0(
            seqnames(junc), ":",
            start(junc), "-", end(junc), "/", strand(junc)
        ),
        juncStart = start(junc) - 1,
        juncEnd = end(junc) + 1
    )

    final <- c()
    isMean <- ("mean" %in% names(mcols(junc)))
    for(i in seq_len(length(junc))) {
        count <- countsd <- 0
        if(isMean) {
            count <- mcols(junc)$mean[i]
            countsd <- mcols(junc)$sd[i]
        } else {
            count <- mcols(junc)$count[i]
        }
        if(count >= junctionThreshold) {
            leftY <- mcols(junc)$leftCoordHeight[i]
            rightY <- mcols(junc)$rightCoordHeight[i]
            leftX <- df$juncStart[i]
            rightX <- df$juncEnd[i]
            juncName <- df$juncName[i]
            
            outdf <- data.frame(
                x = seq(leftX, rightX, length.out = 90),
                y = seq(leftY, rightY, length.out = 90)
            )
            outdf$yarc <- outdf$y + sinpi(seq(0,1,length.out = 90)) * arcHeight
            outdf$coord <- juncName

            if(isMean) {
                outdf$value <- paste0(
                    round(100 * count, 1), "+/-", 
                    round(100 * countsd, 1), " %"
                )
                outdf$info <- paste(
                    paste0("Junction: ", juncName),
                    paste0("PSI: ", outdf$value),
                    sep = "\n"
                )
            } else {
                outdf$value <- count            
                outdf$info <- paste(
                    paste0("Junction: ", juncName),
                    paste0("Depth: ", count),
                    sep = "\n"
                )
            }
            final <- rbind(final, outdf)
        }
    }
    return(final)
}
