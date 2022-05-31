#' @describeIn ASEFilter-class Constructs a ASEFilter object
#' @export
ASEFilter <- function(
        filterClass = c("Data", "Annotation"),
        filterType = c(
            "Depth", "Participation", "Consistency",
            "Modality", "Protein_Coding", "NMD", "TSL", "Terminus", 
            "ExclusiveMXE"
        ),
        pcTRUE = 100, minimum = 20, maximum = 1, minDepth = 5,
        condition = "", minCond = -1,
        EventTypes = c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")
) {
    nfv <- new("ASEFilter",
        filterClass = filterClass, filterType = filterType, pcTRUE = pcTRUE,
        minimum = minimum, maximum = maximum, minDepth = minDepth,
        condition = condition, minCond = minCond,
        EventTypes = EventTypes
    )
    nfv
}

setMethod("initialize", "ASEFilter", function(.Object,
        filterClass = c("Data", "Annotation"),
        filterType = c(
            "Depth", "Participation", "Consistency",
            "Modality", "Protein_Coding", "NMD", "TSL", "Terminus", 
            "ExclusiveMXE"
        ),
        pcTRUE = 100, minimum = 20, maximum = 1, minDepth = 5,
        condition = "", minCond = -1,
        EventTypes = c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")
) {
    .Object <- callNextMethod()

    # filterClass <- match.arg(filterClass)
    # filterType <- match.arg(filterType)

    data_methods <- c("Depth", "Participation", "Consistency")
    annotation_methods <- c("Modality", "Protein_Coding", "NMD", "TSL",
        "Terminus", "ExclusiveMXE")
    
    if(filterClass %in% c("Data", "Annotation")) {
        filterClass <- filterClass[1]
    } else {
        filterClass <- "(none)"
    }
    filterType <- filterType[1]
    
    if (filterClass == "Data") {
        if (!(filterType %in% data_methods))
        .log(paste("filterClass must be a recognised Data method",
            paste(data_methods, collapse = ", ")))
    } else if (filterClass == "Annotation") {
        if (!(filterType %in% annotation_methods))
            .log(paste("filterClass must be a recognised Annotation method",
                paste(annotation_methods, collapse = ", ")))
    } else {
        filterType <- "(none)"
    }

    pcTRUE <- as.numeric(pcTRUE)
    pcTRUE <- min(100, max(0, pcTRUE))
    minimum <- as.numeric(minimum)
    minimum <- max(0, minimum)
    maximum <- as.numeric(maximum)
    maximum <- max(0, maximum)
    minDepth <- as.numeric(minDepth)
    minDepth <- max(0, minDepth)
    minCond <- as.numeric(minCond)

    # Only 1 condition allowed
    if (length(condition) != 1) .log("`condition` should be a character scalar")

    et_args <- c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")
    EventTypes <- EventTypes[EventTypes %in% et_args]
    if (length(EventTypes) == 0) EventTypes <- et_args

    .Object@filterClass <- filterClass
    .Object@filterType <- filterType
    .Object@pcTRUE <- pcTRUE
    .Object@minimum <- minimum
    .Object@maximum <- maximum
    .Object@minDepth <- minDepth
    .Object@condition <- condition
    .Object@minCond <- minCond
    .Object@EventTypes <- EventTypes

    .Object
})

setMethod("show", "ASEFilter", function(object) {
    .cat_ASEFilter_common(object)
    .cat_ASEFilter_filterSpecific(object)
})

# Describe class, type, and conditions filter will be applied across
.cat_ASEFilter_common <- function(object) {
    .nxtcat(paste0("ASEFilter Class: ", .colourise("%s\t", "red")),
        object@filterClass)
    .nxtcat(paste0("Type: ", .colourise("%s\n", "purple")),
        object@filterType)
    .nxtcat(paste0("EventTypes: ", .colourise("%s\n", "green")),
        paste(object@EventTypes, collapse = " "))
    if (object@filterClass == "Data") {
        cat("Filter must be passed in at least ")
        if (object@condition == "") {
            .nxtcat(paste0(.colourise("%.1f", "purple"),
                " percent of all samples\n"), object@pcTRUE)
        } else {
            if (object@minCond == -1) {
            .nxtcat(paste0(.colourise("%.1f", "purple"),
                " percent of all categories of ",
                .colourise("%s\n", "purple")), object@pcTRUE, object@condition)
            } else {
            .nxtcat(paste0(.colourise("%.1f", "purple"),
                " percent of ", .colourise("%i", "green"),
                " categories of ", .colourise("%s\n", "purple")),
                object@pcTRUE, as.integer(object@minCond), object@condition)
            }
        }
    }
}

# Describe filter-specific functions
.cat_ASEFilter_filterSpecific <- function(object) {
    if (object@filterType == "Depth") {
        .nxtcat(paste0("Minimum Event Depth: ",
            .colourise("%i\n", "red")), as.integer(object@minimum))
        .cat_filter_info("minDepth")
    } else if (object@filterType == "Participation") {
        .nxtcat(paste0("Minimum Coverage (Participation) Percentage: ",
            .colourise("%.1f\n", "red")), object@minimum)
        .cat_filter_info("Coverage (Participation) minimum", object@EventTypes)
        .nxtcat(paste0("Event Depth below ", .colourise("%i", "purple"),
            " are ignored\n"), as.integer(object@minDepth))
        .cat_filter_info("minDepth")
    } else if (object@filterType == "Consistency") {
        .nxtcat(paste0("Sub-event minimum proportion: ",
            .colourise("%.1f", "red"), " of total isoform splice count\n"),
            2^(-object@maximum))
        .cat_filter_info("Consistency maximum", object@EventTypes)
        .nxtcat(paste0("Event Depth below ", .colourise("%i", "purple"),
            " are ignored\n"), as.integer(object@minDepth))
        .cat_filter_info("minDepth")
    } else if (object@filterType == "Modality") {
        cat("Events of the following modality are retained")
    } else if (object@filterType == "Protein_Coding") {
        cat("Events of which neither isoform encodes protein are removed")
    } else if (object@filterType == "NMD") {
        cat("Events of which neither isoform are NMD substrates are removed")
    } else if (object@filterType == "TSL") {
        .nxtcat(paste0("Transcript Support Level: ",
            .colourise("%i", "purple"), " or higher\n"),
            as.integer(object@minimum))
        cat("Events with both isoforms belonging to ")
        cat("lower-ranking TSLs are removed")
    } else if (object@filterType == "Terminus") {
        cat("For alternative first / last exons, events overlapping with non-")
        cat("first/last introns are removed")
    } else if (object@filterType == "ExclusiveMXE") {
        cat("MXE events with overlapping casette exons are removed")
    }
}

# Describe more info re specific functions
.cat_filter_info <- function(mode, EventTypes = "") {
    if (mode == "minDepth") {
        cat("Event Depth refers to number of aligned splice reads plus ")
        cat("effective depth of their introns\n")
    } else if (mode == "Coverage (Participation) minimum") {
        if (any(EventTypes %in% c("IR", "RI"))) {
            cat("In retained introns, participation refers to the proportion ")
            cat("of the measured intron that is covered by at least ")
            cat("1 alignment\n")
        }
        if (any(EventTypes %in% c("MXE", "SE", "ALE", "AFE", "A3SS", "A5SS"))) {
            cat("In splice events, participation refers to the proportion of ")
            cat("junction reads that belong to either the included or excluded")
            cat(" isoforms of the given event\n")
        }
    } else if (mode == "Consistency maximum") {
        if (any(EventTypes %in% c("IR", "RI"))) {
            cat("In retained introns, sub-events refer to 5'- and 3'-")
            cat("exon-intron spanning alignments\n")
        }
        if (any(EventTypes %in% c("MXE", "SE"))) {
            cat("In mutually exclusive exons and skipped exons, ")
            cat("sub-events refer to the upstream and downstream splice ")
            cat("alignments\n")
        }
    }
}
