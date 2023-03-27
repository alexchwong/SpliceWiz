



.gCD_plotRef <- function(
    DTplotlist, args,    
    reverseGenomeCoords = FALSE,
    add_information = TRUE,
    interactive = FALSE
) {
    view_chr <- args[["view_chr"]]
    view_start <- args[["view_start"]]
    view_end <- args[["view_end"]]

    group.DT <- copy(DTplotlist$group.DT)

    # use one paste function
    if(add_information) {
        group.DT[, c("left_extra", "right_extra") := list("","")]
        if(reverseGenomeCoords) {
            group.DT[get("strand") == "-", c("right_extra") := " -->"]
            group.DT[get("strand") == "+", c("left_extra") := "<-- "]
        } else {
            group.DT[get("strand") == "+", c("right_extra") := " -->"]
            group.DT[get("strand") == "-", c("left_extra") := "<-- "]
        }
        group.DT[, c("display_name") := paste0(
            get("left_extra"), 
            get("group_name"),
            " (", get("strand"), ") ",
            get("group_biotype"),
            get("right_extra")
        )]
        group.DT[, c("left_extra", "right_extra") := list(NULL,NULL)]

        group.DT[, c("disp_x") := 0.5 * (get("start") + get("end"))]
        group.DT[get("start") < view_start & get("end") > view_start,
            c("disp_x") := 0.5 * (view_start + get("end"))]
        group.DT[get("end") > view_end & get("start") < view_end,
            c("disp_x") := 0.5 * (get("start") + view_end)]
        group.DT[get("start") < view_start & get("end") > view_end,
            c("disp_x") := 0.5 * (view_start + view_end)]    
    }

    reduced <- copy(DTplotlist$reduced.DT)
    reduced <- reduced[!is.na(reduced$plot_level)]
    condense_this <- DTplotlist$condense_this
    exonRanges <- DTplotlist$exonRanges

    # Hover Text annotation
    reduced[, c("Information") := ""]
    reduced[get("type") %in% c("exon", "CDS"), c("Information") := paste(
        paste(get("transcript_id"), "exon", get("aux_id")),
        paste0("(", get("feature_id"), ")"),
        paste0(get("seqnames"), ":", get("start"), "-", get("end"), "/", 
            get("strand")),
        sep = "\n"
    )]
    reduced[get("type") == "intron", c("Information") := paste(
        get("feature_id"), 
        paste0(get("seqnames"), ":", get("start"), "-", get("end"), "/", 
            get("strand")),
        sep = "\n"
    )]
    
    reduced <- as.data.frame(reduced)
    p <- suppressWarnings(ggplot(reduced, aes(text = get("Information"))))
    
    if (nrow(subset(reduced, type = "intron")) > 0) {
        reducedIntrons <- reduced[reduced$type == "intron", ]
        reducedIntronsExpanded <- c()
        for(i in seq_len(nrow(reducedIntrons))) {
            reducedIntronsExpanded <- rbind(reducedIntronsExpanded, data.frame(
                start = seq(reducedIntrons$start[i], reducedIntrons$end[i],
                    length.out = 10),
                end = seq(reducedIntrons$start[i], reducedIntrons$end[i],
                    length.out = 10),
                plot_level = reducedIntrons$plot_level[i],
                highlight = reducedIntrons$highlight[i],
                Information = reducedIntrons$Information[i]              
            ))
        }
        p <- p + geom_line(data = reducedIntronsExpanded,
            aes(x = get("start"), y = get("plot_level"),
            color = get("highlight"), group = get("Information")))
    }
    if (nrow(reduced[reduced$type != "intron", ]) > 0) {
        p <- p + geom_rect(data = reduced[reduced$type != "intron", ],
            aes(xmin = get("start"), xmax = get("end"),
                ymin = get("plot_level") - 0.1 -
                    ifelse(get("type") %in%
                        c("CDS", "start_codon", "stop_codon"), 0.1, 0),
                ymax = get("plot_level") + 0.1 +
                    ifelse(get("type") %in%
                        c("CDS", "start_codon", "stop_codon"), 0.1, 0),
                fill = get("highlight")
            )
        )
    }
    
    col_highlights <- sort(unique(DTplotlist$reduced.DT$highlight))
    col_highlights <- sub("0", "black", col_highlights)
    col_highlights <- sub("1", "blue", col_highlights)
    col_highlights <- sub("2", "red", col_highlights)
    p <- p + scale_color_manual(values = col_highlights) +
        scale_fill_manual(values = col_highlights)    
    
    p <- p + theme_white_legend_plot_track +
        theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
            legend.title = element_blank())
    
    anno <- NULL
    if(add_information) {
        if (condense_this == TRUE) {
            anno <- list(
                x = group.DT$disp_x,
                y = group.DT$plot_level - 0.5 + 0.3 * 
                    runif(rep(1, nrow(group.DT))),
                text = group.DT$display_name,
                xref = "x", yref = "y", showarrow = FALSE)
        } else {
            anno <- list(
                x = group.DT$disp_x,
                y = group.DT$plot_level - 0.4,
                text = group.DT$display_name,
                xref = "x", yref = "y", showarrow = FALSE)
        }
    }

    if (nrow(group.DT) == 0) {
        max_plot_level <- 1
    } else {
        max_plot_level <- max(group.DT$plot_level)
    }

    if(!interactive) {
        out_p <- p
        if(!is.null(anno)) {
            out_p <- out_p + geom_text(
                data = data.frame(
                    x = anno[["x"]], y = anno[["y"]],
                    Information = anno[["text"]]
                ),
                    aes(x = get("x"), y = get("y"), label = get("Information"))
                )
        }
            
        out_p <- out_p + theme(legend.position = "none") +
            labs(x = paste("Chromosome", view_chr))
        ref_ymin <- min(layer_scales(out_p)$y$range$range)
        ref_ymax <- max(layer_scales(out_p)$y$range$range)
        if(!reverseGenomeCoords) {
            plotViewStart <- view_start
            plotViewEnd <- view_end
        } else {
            plotViewStart <- view_end
            plotViewEnd <- view_start     
        }
        out_p <- out_p + 
            scale_x_continuous(labels = label_number(scale_cut = cut_si(""))) +
            coord_cartesian(
                xlim = c(plotViewStart, plotViewEnd),
                ylim = c(ref_ymin - 1, ref_ymax + 1),
                expand = FALSE
            )
    } else {
        if(!reverseGenomeCoords) {
            out_p <- ggplotly(p, tooltip = "text") %>%
                layout(
                    annotations = anno, dragmode = "pan",
                    xaxis = list(range = c(view_start, view_end),
                        title = paste("Chromosome/Scaffold", view_chr)),
                    yaxis = list(range = c(0, 1 + max_plot_level),
                        fixedrange = TRUE)
                ) 
        } else {
            out_p <- ggplotly(p, tooltip = "text") %>%
                layout(
                    annotations = anno, dragmode = "pan",
                    xaxis = list(
                        range = c(view_end, view_start),
                        title = paste("Chromosome/Scaffold", view_chr)
                    ),
                    yaxis = list(range = c(0, 1 + max_plot_level),
                        fixedrange = TRUE)
                )         
        }
        for (i in seq_len(length(out_p$pl$x$data))) {
            out_p$pl$x$data[[i]]$showlegend <- FALSE
        }    
    }
    
    return(out_p)
}

################################################################################

# Get normalized coverage from obj
.gCD_getCoverage <- function(obj, args, normalize = FALSE) {
    view_start <- args[["view_start"]]
    view_end <- args[["view_end"]]
    view_strand <- args[["view_strand"]]
    Event <- args[["Event"]]
    
    outList <- list()
    sampleList <- args[["sampleList"]]
    
    if(view_strand == "+") {
        cov <- obj@covData$pos
    } else if(view_strand == "-") {
        cov <- obj@covData$neg    
    } else {
        cov <- obj@covData$uns    
    }
    
    for(i in seq_len(length(sampleList))) {
        track <- names(sampleList)[i]
        samples <- sampleList[[i]]
        
        matList <- list()
        for(s in samples) {
            view <- IRanges::Views(cov[[s]], view_start, view_end)
            matList[[s]] <- as.matrix(view[[1]])
        }
        mat <- do.call(cbind, matList)
        colnames(mat) <- samples
        
        # Normalize coverage
        if(normalize) {
            normFactor <- unlist(obj@normData$norms[Event,samples])
            mat <- t(t(mat) / normFactor)
        }
        
        outList[[track]] <- cbind(seq(view_start, view_end), as.data.frame(mat))
        colnames(outList[[track]])[1] <- "x"
        
        if(ncol(outList[[track]]) == 2) {
            colnames(outList[[track]])[2] <- "depth"
        }
    }
    
    return(outList)
}

.gCD_binCoverage <- function(cov, divFactor, junc) {
    out <- list()
    outnames <- names(cov)
    
    juncSpares <- c()
    for(i in seq_len(length(junc))) {
        juncSpares <- c(
            juncSpares,
            start(junc[[i]]) - 1,
            end(junc[[i]]) + 1
        )
    }
    juncSpares <- sort(unique(juncSpares))
    
    for(i in seq_len(length(cov))) {
        out[[outnames[i]]] <- .gCD_bin_df(
            cov[[outnames[i]]], 
            divFactor, juncSpares
        )
    }
    return(out)
}

.gCD_bin_df <- function(df, binwidth = 3, spareCoords = c()) {
    DT <- as.data.table(df)
    brks <- seq(1, 
        nrow(DT) + 1, 
        length.out = (nrow(DT) + 1) / binwidth
    )
    if(length(spareCoords) > 0) {
        brks <- c(brks, which(DT$x %in% spareCoords))
    }
    brks <- sort(unique(brks))
    
    bin <- NULL
    DT[, c("bin") := findInterval(seq_len(nrow(DT)), brks)]
    DT2 <- DT[, lapply(.SD, mean, na.rm = TRUE), by = "bin"]
    DT2[, c("bin") := NULL]
    return(as.data.frame(DT2))
}

.gCD_covStats <- function(cov, ribbon_mode) {
    out <- list()
    outnames <- names(cov)
    
    for(i in seq_len(length(cov))) {
        covmat <- as.matrix(cov[[outnames[i]]][,-1])
        out[[outnames[i]]] <- data.frame(
            x = cov[[outnames[i]]]$x,
            mean = rowMeans(covmat)
        )
        n <- ncol(covmat)
        if(ribbon_mode == "sd") {
            out[[outnames[i]]]$var <- rowSds(covmat)
        } else if(ribbon_mode == "sem") {
            out[[outnames[i]]]$var <- rowSds(covmat) / sqrt(n)
        } else if(ribbon_mode == "ci") {
            conf.int <- 0.95
            out[[outnames[i]]]$var <- qt((1 + conf.int) / 2, df = n - 1) * 
                rowSds(covmat) / sqrt(n)
        }
    }
    return(out)
}

.gCD_ttest <- function(cov1, cov2) {
    coords <- sort(intersect(cov1$x, cov2$x))
    mat <- cbind(as.matrix(
        cov1[cov1$x %in% coords, -1],
        cov2[cov2$x %in% coords, -1]
    ))
    fac <- factor(rep(c("1", "2"), each = c(ncol(cov1) - 1, ncol(cov2) - 1)))
    
    t_test <- genefilter::rowttests(mat, fac)
    
    ret <- data.frame(
        x = coords,
        t_stat = -log10(t_test$p.value)
    )
    ret$t_stat[!is.finite(ret$t_stat)] <- 0
    return(ret)
}

.gCD_getPSI <- function(obj, args, norm_cov) {
    view_chr <- args[["view_chr"]]
    view_start <- args[["view_start"]]
    view_end <- args[["view_end"]]
    view_strand_jn <- args[["view_strand_jn"]]

    outList <- list()
    sampleList <- args[["sampleList"]]

    gr_select <- GRanges(
        view_chr, 
        IRanges(view_start, view_end), 
        view_strand_jn
    )
    OL <- findOverlaps(obj@juncData$junc_gr, gr_select)
        
    if(length(unique(from(OL))) == 0) return(NULL)

    hits <- sort(unique(from(OL)))
    junc_gr <- obj@juncData$junc_gr[hits]
    leftCoords <- start(junc_gr) - 1
    rightCoords <- end(junc_gr) + 1
    lC_inrange <- leftCoords >= view_start & leftCoords <= view_end
    rC_inrange <- rightCoords >= view_start & rightCoords <= view_end

    # normalize by starting coordinates:

    for(i in seq_len(length(sampleList))) {
        track <- names(sampleList)[i]
        samples <- sampleList[[i]]
        
        junc <- obj@juncData$junc_PSI[hits, samples, drop = FALSE]
        
        # Returned data
        outList[[track]] <- .grDT(cbind(
            as.data.frame(junc_gr),
            data.frame(
                mean = rowMeans(junc),
                sd = rowSds(junc),
                leftCoordHeight = 1,
                rightCoordHeight = 1
            )
        ), keep.extra.columns = TRUE)
        
        cov <- norm_cov[[track]]
        
        if(sum(lC_inrange) > 0) {
            covL <- cov[leftCoords[lC_inrange] - cov$x[1] + 1, -1]
            mcols(outList[[track]])$leftCoordHeight[lC_inrange] <- rowMeans(covL)
        }
            
        if(sum(rC_inrange) > 0) {
            covR <- cov[rightCoords[rC_inrange] - cov$x[1] + 1, -1]
            mcols(outList[[track]])$rightCoordHeight[rC_inrange] <- rowMeans(covR)
        }
    }
    return(outList)
}

.gCD_getJuncRaw <- function(obj, args, raw_cov) {
    view_chr <- args[["view_chr"]]
    view_start <- args[["view_start"]]
    view_end <- args[["view_end"]]
    view_strand_jn <- args[["view_strand_jn"]]

    outList <- list()
    sampleList <- args[["sampleList"]]

    gr_select <- GRanges(
        view_chr, 
        IRanges(view_start, view_end), 
        view_strand_jn
    )
    OL <- findOverlaps(obj@juncData$junc_gr, gr_select)
        
    if(length(unique(from(OL))) == 0) return(NULL)

    hits <- sort(unique(from(OL)))
    junc_gr <- obj@juncData$junc_gr[hits]
    leftCoords <- start(junc_gr) - 1
    rightCoords <- end(junc_gr) + 1
    lC_inrange <- leftCoords >= view_start & leftCoords <= view_end
    rC_inrange <- rightCoords >= view_start & rightCoords <= view_end

    # normalize by starting coordinates:

    for(i in seq_len(length(sampleList))) {
        track <- names(sampleList)[i]
        samples <- sampleList[[i]]
        
        if(view_strand_jn == "*") {
            junc <- obj@juncData$junc_counts_uns[hits, samples, drop = FALSE]
        } else {
            junc <- obj@juncData$junc_counts[hits, samples, drop = FALSE]            
        }
        
        cov <- raw_cov[[track]]

        # Returned data
        outList[[track]] <- .grDT(cbind(
            as.data.frame(junc_gr),
            data.frame(
                count = unname(unlist(junc)),
                leftCoordHeight = max(cov[,2]),
                rightCoordHeight = max(cov[,2])
            )
        ), keep.extra.columns = TRUE)

        if(sum(lC_inrange) > 0) {
            covL <- cov[leftCoords[lC_inrange] - cov$x[1] + 1, -1]
            mcols(outList[[track]])$leftCoordHeight[lC_inrange] <- covL
        }
        if(sum(rC_inrange) > 0) {
            covR <- cov[rightCoords[rC_inrange] - cov$x[1] + 1, -1]
            mcols(outList[[track]])$rightCoordHeight[rC_inrange] <- covR
        }
    }
    return(outList)
}

