# Explicitly define constructor
setMethod("initialize", "NxtSE", 
    function(.Object, ...) {
        .Object@int_elementMetadata <- S4Vectors::DataFrame()
        .Object@int_colData <- S4Vectors::DataFrame()
        
        .Object <- callNextMethod()

        .Object
    })

#' @describeIn NxtSE-class Constructor function for NxtSE; akin to
#'   SummarizedExperiment(...)
NxtSE <- function(...) {
    se <- SummarizedExperiment(...)

    .se_to_nxtse(se)
}

#' @exportMethod coerce
setAs("SummarizedExperiment", "NxtSE", function(from) {
    .se_to_nxtse(from)
})

# Converts a SummarizedExperiment to NxtSE object. Adapted from SingleCellExp
.se_to_nxtse <- function(se) {
    old <- S4_disableValidity()
    if (!isTRUE(old)) {
        S4_disableValidity(TRUE)
        on.exit(S4_disableValidity(old))
    }

    if (!all(c("Included", "Excluded", "Depth", "Coverage", "minDepth") %in%
        names(assays(se)))) {
        .log(paste("Object was not created by makeSE(),",
            "returning SummarizedExperiment object instead"), "warning",
            use_system_time = FALSE)
        return(se)
    }
    if (!all(c("Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc") %in%
            names(S4Vectors::metadata(se)))) {
        .log(paste("Object was not created by makeSE(),",
            "returning SummarizedExperiment object instead"), "warning",
            use_system_time = FALSE)
        return(se)
    }
    out <- new("NxtSE", se)
    up_inc(out) <- se@metadata[["Up_Inc"]]
    down_inc(out) <- se@metadata[["Down_Inc"]]
    up_exc(out) <- se@metadata[["Up_Exc"]]
    down_exc(out) <- se@metadata[["Down_Exc"]]
    covfile(out) <- se@metadata[["cov_file"]]
    sampleQC(out) <- se@metadata[["sampleQC"]]
    ref(out) <- se@metadata[["ref"]]
    return(out)
}

# Check validity of metadata "Inc"
.valid.NxtSE.meta_inc <- function(x) {
    if (!("Up_Inc" %in% names(metadata(x)))) {
        txt <- "Up_Inc is not found in NxtSE object"
        return(txt)
    }
    if (!("Down_Inc" %in% names(metadata(x)))) {
        txt <- "Down_Inc is not found in NxtSE object"
        return(txt)
    }
    rownames_rowData <- rowData(x)$EventName[
        rowData(x)$EventType %in% c("IR", "SE", "MXE")
    ]
    rownames_up_inc <- rownames(up_inc(x))
    if (!identical(rownames_rowData, rownames_up_inc)) {
        txt <- paste("Up_Inc events do not encompass all IR / SE / MXE events",
            "in rowData")
        return(txt)
    }
    rownames_down_inc <- rownames(down_inc(x))
    if (!identical(rownames_rowData, rownames_down_inc)) {
        txt <- paste("Down_Inc events do not encompass all IR / SE / MXE",
            "events in rowData")
        return(txt)
    }
    NULL
}

# Check validity of metadata "Exc"
.valid.NxtSE.meta_exc <- function(x) {
    if (!("Up_Exc" %in% names(metadata(x)))) {
        txt <- "Up_Exc is not found in NxtSE object"
        return(txt)
    }
    if (!("Down_Exc" %in% names(metadata(x)))) {
        txt <- "Down_Exc is not found in NxtSE object"
        return(txt)
    }
    rownames_rowData <- rowData(x)$EventName[
        rowData(x)$EventType %in% c("MXE")
    ]
    rownames_up_exc <- rownames(up_exc(x))
    if (!identical(rownames_rowData, rownames_up_exc)) {
        txt <- paste("Up_exc events do not encompass all IR / SE /",
            "MXE events in rowData")
        return(txt)
    }
    rownames_down_exc <- rownames(down_exc(x))
    if (!identical(rownames_rowData, rownames_down_exc)) {
        txt <- paste("Down_exc events do not encompass all IR / SE /",
            "MXE events in rowData")
        return(txt)
    }
    NULL
}

# Check validity of metadata "COV". Invalidates if some files are not COV
# Will validate if all files are COV, or if all are empty
.valid.NxtSE.meta_cov <- function(x) {
    cov_files <- metadata(x)[["cov_file"]]
    if (!all(file.exists(cov_files) | cov_files == "")) {
        txt <- paste(
            "Some coverage files are not found:",
            paste(cov_files[!file.exists(cov_files)], collapse = " ")
        )
        return(txt)
    }
    if (!isCOV(cov_files)) {
        txt <- "Some coverage files are not validated"
        return(txt)
    }
    if (length(cov_files) != ncol(x)) {
        txt <- "cov_files must be of same length as number of samples"
        return(txt)
    }
    if (!identical(names(cov_files), colnames(x))) {
        txt <- paste("names of cov_files vector must be identical to",
            "colnames of NxtSE object")
        return(txt)
    }
    NULL
}

# Check validity of QC statistics
.valid.NxtSE.meta_QC <- function(x) {
    sampleQC <- metadata(x)[["sampleQC"]]
    if (!identical(rownames(sampleQC), colnames(x))) {
        txt <- paste("row names of sampleQC data frame must be identical to",
            "colnames of NxtSE object")
        return(txt)
    }
    NULL
}

# Main validity check @ metadata
.valid.NxtSE.metadata <- function(x)
{
    c(
        .valid.NxtSE.meta_inc(x),
        .valid.NxtSE.meta_exc(x),
        .valid.NxtSE.meta_cov(x),
        .valid.NxtSE.meta_QC(x)
    )
}

.valid.NxtSE <- function(x)
{
    .valid.NxtSE.metadata(x)
}

setValidity2("NxtSE", .valid.NxtSE)

################################################################################
# Convenience functions. Some are adapted from SingleCellExperiment

setMethod("show", "NxtSE",
    function(object)
{
    callNextMethod(object)
})

simplify_NULL_dimnames <- function(dimnames)
{
    if (all(vapply(dimnames, is.null, logical(1), USE.NAMES = FALSE))) {
        return(NULL)
    }
    dimnames
}

set_dimnames <- function(x, value)
{
    if (!identical(dimnames(x), value)) {
dimnames(x) <- value
}
    x
}

assay_withDimnames <- function(x, assay) {
    x_dimnames <- dimnames(x)
    if (is.null(x_dimnames)) {
x_dimnames <- list(NULL, NULL)
}
    a_dimnames <- dimnames(assay)
    x_dimnames[[1]] <- x_dimnames[[1]][x_dimnames[[1]] %in% a_dimnames[[1]]]
    a_dimnames[c(1, 2)] <- x_dimnames
    a_dimnames <- simplify_NULL_dimnames(a_dimnames)
    set_dimnames(assay, a_dimnames)
}

vector_withDimnames <- function(x, vector) {
    x_dimnames <- dimnames(x)[[2]]
    if (is.null(names(vector)) || any(is.na(names(vector)))) {
        names(vector) <- x_dimnames
        return(vector)
    }
    if (is.null(x_dimnames) || any(is.na(x_dimnames))) {
        .log("Some sample names in NxtSE have null values")
    }
    vector[x_dimnames]
}

sampleQC_withDimnames <- function(x, df) {
    x_dimnames <- dimnames(x)[[2]]
    if (is.null(rownames(df)) || any(is.na(rownames(df)))) {
        rownames(df) <- x_dimnames
        return(df)
    }
    if (is.null(x_dimnames) || any(is.na(x_dimnames))) {
        .log("Some sample names in NxtSE have null values")
    }
    df[x_dimnames, , drop = FALSE]
}

#' @describeIn NxtSE-class Gets upstream included events (SE/MXE), or
#'   upstream exon-intron spanning reads (IR)
#' @export
setMethod("up_inc", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    if (withDimnames) {
        assay_withDimnames(x, x@metadata[["Up_Inc"]][, , drop = FALSE])
    } else {
        x@metadata[["Up_Inc"]][, , drop = FALSE]
    }
})

################################## GETTERS #####################################

#' @describeIn NxtSE-class Gets downstream included events (SE/MXE), or
#'   downstream exon-intron spanning reads (IR)
#' @export
setMethod("down_inc", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    if (withDimnames) {
        assay_withDimnames(x, x@metadata[["Down_Inc"]][, , drop = FALSE])
    } else {
        x@metadata[["Down_Inc"]][, , drop = FALSE]
    }
})

#' @describeIn NxtSE-class Gets upstream excluded events (MXE only)
#' @export
setMethod("up_exc", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    if (withDimnames) {
        assay_withDimnames(x, x@metadata[["Up_Exc"]][, , drop = FALSE])
    } else {
        x@metadata[["Up_Exc"]][, , drop = FALSE]
    }
})

#' @describeIn NxtSE-class Gets downstream excluded events (MXE only)
#' @export
setMethod("down_exc", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    if (withDimnames) {
        assay_withDimnames(x, x@metadata[["Down_Exc"]][, , drop = FALSE])
    } else {
        x@metadata[["Down_Exc"]][, , drop = FALSE]
    }
})

#' @describeIn NxtSE-class Gets a named vector with
#'   the paths to the corresponding COV files
#' @export
setMethod("covfile", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    if (withDimnames) {
        vector_withDimnames(x, x@metadata[["cov_file"]])
    } else {
        x@metadata[["cov_file"]]
    }
})

#' @describeIn NxtSE-class Gets a data frame with the QC parameters
#'   of the samples
#' @export
setMethod("sampleQC", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    if (withDimnames) {
        sampleQC_withDimnames(x, x@metadata[["sampleQC"]])
    } else {
        x@metadata[["sampleQC"]]
    }
})

#' @describeIn NxtSE-class Retrieves a list of annotation data associated
#'   with this NxtSE object; primarily used in plotCoverage()
#' @export
setMethod("ref", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    x@metadata[["ref"]]
})

#' @describeIn NxtSE-class Converts all DelayedMatrix assays as matrices
#'   (i.e. performs all delayed calculation and loads resulting object
#'   to RAM)
#' @export
setMethod("realize_NxtSE", c("NxtSE"), function(x, withDimnames = TRUE, ...) {
    assay.todo <- c("Included", "Excluded", "Depth", "Coverage",
        "minDepth")
    for (assayname in assay.todo) {
        assay(x, assayname) <- as.matrix(assay(x, assayname))
    }
    up_inc(x) <- as.matrix(up_inc(x))
    down_inc(x) <- as.matrix(down_inc(x))
    up_exc(x) <- as.matrix(up_exc(x))
    down_exc(x) <- as.matrix(down_exc(x))
    return(x)
})

################################## SETTERS #####################################

#' @describeIn NxtSE-class Sets upstream included events (SE/MXE), or
#'   upstream exon-intron spanning reads (IR)
#' @export
setReplaceMethod("up_inc", c("NxtSE"), function(x, withDimnames = TRUE, value)
{
    if (!is.matrix(value) & !is(value, "DelayedMatrix")) {
        .log("replacement 'up_inc' value must be a matrix or DelayedMatrix")
    }
    if (withDimnames) {
        value <- assay_withDimnames(x, value)
    }
    x@metadata[["Up_Inc"]] <- value
    x
})

#' @describeIn NxtSE-class Sets downstream included events (SE/MXE), or
#'   downstream exon-intron spanning reads (IR)
#' @export
setReplaceMethod("down_inc", c("NxtSE"), function(x, withDimnames = TRUE, value)
{
    if (!is.matrix(value) & !is(value, "DelayedMatrix")) {
        .log("replacement 'down_inc' value must be a matrix or DelayedMatrix")
    }
    if (withDimnames) {
        value <- assay_withDimnames(x, value)
    }
    x@metadata[["Down_Inc"]] <- value
    x
})

#' @describeIn NxtSE-class Sets upstream excluded events (MXE only)
#' @export
setReplaceMethod("up_exc", c("NxtSE"), function(x, withDimnames = TRUE, value)
{
    if (!is.matrix(value) & !is(value, "DelayedMatrix")) {
        .log("replacement 'up_exc' value must be a matrix or DelayedMatrix")
    }
    if (withDimnames) {
        value <- assay_withDimnames(x, value)
    }
    x@metadata[["Up_Exc"]] <- value
    x
})

#' @describeIn NxtSE-class Sets downstream excluded events (MXE only)
#' @export
setReplaceMethod("down_exc", c("NxtSE"), function(x, withDimnames = TRUE, value)
{
    if (!is.matrix(value) & !is(value, "DelayedMatrix")) {
        .log("replacement 'down_exc' value must be a matrix or DelayedMatrix")
    }
    if (withDimnames) {
        value <- assay_withDimnames(x, value)
    }
    x@metadata[["Down_Exc"]] <- value
    x
})

#' @describeIn NxtSE-class Sets the paths to the corresponding COV files
#' @export
setReplaceMethod("covfile", c("NxtSE"), function(x, withDimnames = TRUE, value)
{
    if (withDimnames) {
        value <- vector_withDimnames(x, value)
    }
    x@metadata[["cov_file"]] <- value
    x
})

#' @describeIn NxtSE-class Sets the values in the data frame containing
#'   sample QC
#' @export
setReplaceMethod("sampleQC", c("NxtSE"), function(x, withDimnames = TRUE, value)
{
    if (withDimnames) {
        value <- sampleQC_withDimnames(x, value)
    }
    x@metadata[["sampleQC"]] <- value
    x
})

setReplaceMethod("ref", c("NxtSE"), function(x, value)
{
    x@metadata[["ref"]] <- value
    x
})

################################ SUBSETTERS ####################################

#' @describeIn NxtSE-class Subsets a NxtSE object
#' @export
setMethod("[", c("NxtSE", "ANY", "ANY"), function(x, i, j, ...) {
    old <- S4_disableValidity()
    if (!isTRUE(old)) {
        S4_disableValidity(TRUE)
        on.exit(S4_disableValidity(old))
    }
    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
            i <- SE_charbound(i, rownames(x), fmt)
        }
        ii <- as.vector(i)
    }
    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
            j <- SE_charbound(j, colnames(x), fmt)
        }
        jj <- as.vector(j)
    }
    if (!missing(i) && !missing(j)) {
        events <- rownames(x)[ii]
        events_Inc <- events[events %in% rownames(metadata(x)[["Up_Inc"]])]
        events_Exc <- events[events %in% rownames(metadata(x)[["Up_Exc"]])]
        samples <- colnames(x)[jj]

        up_inc(x, FALSE) <-
            up_inc(x, FALSE)[events_Inc, samples, drop = FALSE]
        down_inc(x, FALSE) <-
            down_inc(x, FALSE)[events_Inc, samples, drop = FALSE]
        up_exc(x, FALSE) <-
            up_exc(x, FALSE)[events_Exc, samples, drop = FALSE]
        down_exc(x, FALSE) <-
            down_exc(x, FALSE)[events_Exc, samples, drop = FALSE]
        covfile(x, FALSE) <- covfile(x, FALSE)[samples]
        sampleQC(x, FALSE) <- sampleQC(x, FALSE)[samples, , drop = FALSE]
    } else if (!missing(i)) {
        events <- rownames(x)[ii]
        events_Inc <- events[events %in% rownames(metadata(x)[["Up_Inc"]])]
        events_Exc <- events[events %in% rownames(metadata(x)[["Up_Exc"]])]

        up_inc(x, FALSE) <-
            up_inc(x, FALSE)[events_Inc, , drop = FALSE]
        down_inc(x, FALSE) <-
            down_inc(x, FALSE)[events_Inc, , drop = FALSE]
        up_exc(x, FALSE) <-
            up_exc(x, FALSE)[events_Exc, , drop = FALSE]
        down_exc(x, FALSE) <-
            down_exc(x, FALSE)[events_Exc, , drop = FALSE]
    } else if (!missing(j)) {
        samples <- colnames(x)[jj]

        up_inc(x, FALSE) <- up_inc(x, FALSE)[, samples, drop = FALSE]
        down_inc(x, FALSE) <- down_inc(x, FALSE)[, samples, drop = FALSE]
        up_exc(x, FALSE) <- up_exc(x, FALSE)[, samples, drop = FALSE]
        down_exc(x, FALSE) <- down_exc(x, FALSE)[, samples, drop = FALSE]
        covfile(x, FALSE) <- covfile(x, FALSE)[samples]
        sampleQC(x, FALSE) <- sampleQC(x, FALSE)[samples, , drop = FALSE]
    }
    callNextMethod()

})

#' @describeIn NxtSE-class Sets a subsetted NxtSE object
#' @export
setMethod("[<-", c("NxtSE", "ANY", "ANY", "NxtSE"),
        function(x, i, j, ..., value) {
    if (missing(i) && missing(j)) {
        return(value)
    }
    if (!missing(i) && !missing(j)) {
        events <- rownames(x)[i]
        events_Inc <- events[events %in% rownames(metadata(x)[["Up_Inc"]])]
        events_Exc <- events[events %in% rownames(metadata(x)[["Up_Exc"]])]
        samples <- colnames(x)[j]

        up_inc(x)[events_Inc, samples] <- up_inc(value)
        down_inc(x)[events_Inc, samples] <- down_inc(value)
        up_exc(x)[events_Exc, samples] <- up_exc(value)
        down_exc(x)[events_Exc, samples] <- down_exc(value)
        covfile(x)[samples] <- covfile(value)
        sampleQC(x)[samples, ] <- sampleQC(value)
    } else if (!missing(i)) {
        events <- rownames(x)[i]
        events_Inc <- events[events %in% rownames(metadata(x)[["Up_Inc"]])]
        events_Exc <- events[events %in% rownames(metadata(x)[["Up_Exc"]])]

        up_inc(x)[events_Inc, ] <- up_inc(value)
        down_inc(x)[events_Inc, ] <- down_inc(value)
        up_exc(x)[events_Exc, ] <- up_exc(value)
        down_exc(x)[events_Exc, ] <- down_exc(value)
    } else if (!missing(j)) {
        samples <- colnames(x)[j]
        up_inc(x)[, samples] <- up_inc(value)
        down_inc(x)[, samples] <- down_inc(value)
        up_exc(x)[, samples] <- up_exc(value)
        down_exc(x)[, samples] <- down_exc(value)
        covfile(x)[samples] <- covfile(value)
        sampleQC(x)[samples, ] <- sampleQC(value)
    }

    callNextMethod()
})

################################# COMBINERS ####################################


#' @describeIn NxtSE-class Combines two NxtSE objects (by samples - columns)
#' @export
setMethod("cbind", "NxtSE", function(..., deparse.level = 1) {
    old <- S4_disableValidity()
    if (!isTRUE(old)) {
        S4_disableValidity(TRUE)
        on.exit(S4_disableValidity(old))
    }
    out <- callNextMethod()

    metadata <- list()
    args <- list(...)
    tryCatch({
        metadata$Up_Inc <- do.call(cbind, lapply(args, up_inc))
    }, error = function(err) {
        stop(
            "failed to combine 'Up_Inc' in 'cbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$Down_Inc <- do.call(cbind, lapply(args, down_inc))
    }, error = function(err) {
        stop(
            "failed to combine 'Down_Inc' in 'cbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$Up_Exc <- do.call(cbind, lapply(args, up_exc))
    }, error = function(err) {
        stop(
            "failed to combine 'Up_Exc' in 'cbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$Down_Exc <- do.call(cbind, lapply(args, down_exc))
    }, error = function(err) {
        stop(
            "failed to combine 'Down_Exc' in 'cbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$cov_file <- do.call(c, lapply(args, covfile))
    }, error = function(err) {
        stop(
            "failed to combine 'cov_file' in 'cbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$sampleQC <- do.call(rbind, lapply(args, sampleQC))
    }, error = function(err) {
        stop(
            "failed to combine 'sampleQC' in 'cbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$ref <- ref(args[[1]])
    }, error = function(err) {
        stop(
            "failed to combine 'cov_file' in 'cbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    BG_replaceSlots(out, metadata = metadata, check = FALSE)
})

#' @describeIn NxtSE-class Combines two NxtSE objects (by AS/IR events - rows)
#' @export
setMethod("rbind", "NxtSE", function(..., deparse.level = 1) {
    old <- S4_disableValidity()
    if (!isTRUE(old)) {
        S4_disableValidity(TRUE)
        on.exit(S4_disableValidity(old))
    }
    out <- callNextMethod()

    metadata <- list()
    args <- list(...)
    tryCatch({
        metadata$Up_Inc <- do.call(rbind, lapply(args, up_inc))
    }, error = function(err) {
        stop(
            "failed to combine 'Up_Inc' in 'rbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$Down_Inc <- do.call(rbind, lapply(args, down_inc))
    }, error = function(err) {
        stop(
            "failed to combine 'Down_Inc' in 'rbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$Up_Exc <- do.call(rbind, lapply(args, up_exc))
    }, error = function(err) {
        stop(
            "failed to combine 'Up_Exc' in 'rbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$Down_Exc <- do.call(rbind, lapply(args, down_exc))
    }, error = function(err) {
        stop(
            "failed to combine 'Down_Exc' in 'rbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$cov_file <- covfile(args[[1]])
    }, error = function(err) {
        stop(
            "failed to combine 'cov_file' in 'rbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$sampleQC <- sampleQC(args[[1]])
    }, error = function(err) {
        stop(
            "failed to combine 'sampleQC' in 'rbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    tryCatch({
        metadata$ref <- ref(args[[1]])
    }, error = function(err) {
        stop(
            "failed to combine 'ref' in 'rbind(<",
            class(args[[1]]), ">)':\n  ",
            conditionMessage(err))
    })
    BG_replaceSlots(out, metadata = metadata, check = FALSE)
})
