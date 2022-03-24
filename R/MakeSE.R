#' Constructs a NxtSE object from the collated data
#'
#' Creates a \linkS4class{NxtSE} object from the data
#' from IRFinder output collated using [collateData]. This object is used
#' for downstream differential analysis of IR and alternative splicing events
#' using [ASE-methods] as well as visualisation using [plotCoverage]
#'
#' @details
#' `makeSE` retrieves the generic SummarizedExperiment structure saved by
#' [collateData], and initialises a \linkS4class{NxtSE} object. It references
#' the required on-disk assay data using DelayedArrays, thereby utilising
#' 'on-disk' memory to conserve memory usage.
#'
#' For extremely large datasets, loading the entire data into memory may consume
#' too much memory. In such cases, make a subset of the \linkS4class{NxtSE}
#' object (e.g. subset by samples) before loading the data into memory (RAM) 
#' using [realize_NxtSE]
#' 
#' It should be noted that downstream applications of NxtIRF, including
#' [ASE-methods], [plotCoverage], are much faster if the \linkS4class{NxtSE}
#' is realized. It is recommended to realize the \linkS4class{NxtSE} object
#' before extensive usage.
#'
#' If COV files assigned via [collateData] have been moved relative to the
#' `collate_path`, the created \linkS4class{NxtSE} object will not have any
#' linked COV files and [plotCoverage] cannot be used. To reassign these
#' files, a vector of file paths corresponding to all the COV files of the data
#' set can be assigned using `covfile(se) <- vector_of_cov_files`. See
#' example below for details.
#'
#' If `RemoveOverlapping = TRUE`, `makeSE` will try to
#' identify which introns belong to major isoforms, then remove introns of
#' minor introns that overlaps those of major isoforms. Non-overlapping
#' introns are then reassessed iteratively, until all introns are included
#' or excluded in this way. This is important to ensure that overlapping
#' novel IR events are not 'double-counted'.
#'
#' @param collate_path (Required) The output path of [collateData] pointing
#'   to the collated data
#' @param colData (Optional) A data frame containing the sample annotation
#'   information. The first column must contain the sample names.
#'   Omit `colData` to generate a NxtSE object of the whole dataset without
#'   any assigned annotations.
#'   Alternatively, if the names of only a subset of samples are given, then
#'   `makeSE()` will construct the NxtSE object based only on the samples given.
#'   The colData can be set later using `colData()`
#' @param RemoveOverlapping (default = `TRUE`) Whether to filter out overlapping
#'   novel IR events belonging to minor isoforms. See details.
#' @param realize (default = `FALSE`) Whether to load all assay data into
#'   memory. See details
#'
#' @return A \linkS4class{NxtSE} object containing the compiled data in
#' DelayedArrays pointing to the assay data contained in the given
#' `collate_path`
#'
#' @examples
#'
#' # The following code can be used to reproduce the NxtSE object
#' # that can be fetched with SpliceWiz_example_NxtSE()
#'
#' buildRef(
#'     reference_path = file.path(tempdir(), "Reference"),
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#'
#' bams <- SpliceWiz_example_bams()
#' processBAM(bams$path, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "SpliceWiz_Output")
#' )
#'
#' expr <- findSpliceWizOutput(file.path(tempdir(), "SpliceWiz_Output"))
#' collateData(expr,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "Collated_output")
#' )
#'
#' se <- makeSE(collate_path = file.path(tempdir(), "Collated_output"))
#'
#' # "Realize" NxtSE object to load all H5 assays into memory:
#'
#' se <- realize_NxtSE(se)
#'
#' # If COV files have been removed since the last call to collateData()
#' # reassign them to the NxtSE object, for example:
#'
#' covfile_path <- system.file("extdata", package = "SpliceWiz")
#' covfile_df <- findSamples(covfile_path, ".cov")
#'
#' covfile(se) <- covfile_df$path
#'
#' # Check that the produced object is identical to the example NxtSE
#'
#' example_se <- SpliceWiz_example_NxtSE()
#' identical(se, example_se) # should return TRUE
#' @md
#' @export
makeSE <- function(
        collate_path, colData, RemoveOverlapping = TRUE,
        realize = FALSE
) {
    # Includes iterative filtering for IR events with highest mean PSI
        # To annotate IR events of major isoforms

    colData <- .makeSE_validate_args(collate_path, colData)
    colData <- .makeSE_colData_clean(colData)

    collate_path <- normalizePath(collate_path)

    N <- 3
    dash_progress("Loading NxtSE object from file...", N)
    .log("Loading NxtSE object from file...", "message", appendLF = FALSE)

    se <- .makeSE_load_NxtSE(file.path(collate_path, "NxtSE.rds"))

    # Locate relative paths of COV files, or have all-empty if not all are found
    covfiles <- file.path(collate_path, se@metadata[["cov_file"]])
    display_cov_missing_message <- FALSE
    if (all(se@metadata[["cov_file"]] == "") || any(!file.exists(covfiles))) {
        se@metadata[["cov_file"]] <- rep("", ncol(se))
        display_cov_missing_message <- TRUE
    } else {
        se@metadata[["cov_file"]] <- normalizePath(covfiles)
    }

    # Encapsulate as NxtSE object
    se <- as(se, "NxtSE")

    # Subset
    se <- se[, colData$sample]
    if (ncol(colData) > 1) {
        colData_use <- colData[, -1, drop = FALSE]
        rownames(colData_use) <- colData$sample
        colData(se) <- as(colData_use, "DataFrame")
    }

    message("done\n")
    if(display_cov_missing_message)
        .log(paste("Coverage files were not set or not found.",
            "To set coverage files, use `covfile(se) <- filenames`")
        , "message")

    if (realize == TRUE) {
        dash_progress("Realizing NxtSE object...", N)
        .log("Realizing NxtSE object...", "message")
        se <- realize_NxtSE(se)
    }
    
    if (RemoveOverlapping == TRUE) {
        dash_progress("Removing overlapping introns...", N)
        .log("Removing overlapping introns...", "message")
        se <- .makeSE_iterate_IR(se, collate_path)
    }

    # Remove events with NA's (not sure why this occurs)
    Inc_NA <- is.na(rowSums(assay(se, "Included")))
    Exc_NA <- is.na(rowSums(assay(se, "Excluded")))
    
    se <- se[!Inc_NA & !Exc_NA,]

    Up_Inc_NA <- rownames(up_inc(se))[is.na(rowSums(up_inc(se)))]
    Down_Inc_NA <- rownames(down_inc(se))[is.na(rowSums(down_inc(se)))]
    Up_Exc_NA <- rownames(up_exc(se))[is.na(rowSums(up_exc(se)))]
    Down_Exc_NA <- rownames(down_exc(se))[is.na(rowSums(down_exc(se)))]
    
    names_NA <- unique(c(Up_Inc_NA, Down_Inc_NA, Up_Exc_NA, Down_Exc_NA))

    se <- se[!(rownames(se) %in% names_NA),]
    
    return(se)
}



################################################################################

# Checks:
# - whether the given path contains a valid collateData() output
# - whether
.makeSE_validate_args <- function(collate_path, colData) {
    item.todo <- c("rowEvent", "Included", "Excluded", "Depth", "Coverage",
        "minDepth", "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")

    if (!file.exists(file.path(collate_path, "colData.Rds"))) {
        .log(paste("In makeSE():",
            file.path(collate_path, "colData.Rds"),
            "was not found"))
    }
    colData.Rds <- readRDS(file.path(collate_path, "colData.Rds"))
    if (!("df.anno" %in% names(colData.Rds))) {
        .log(paste("In makeSE():",
            file.path(collate_path, "colData.Rds"),
            "must contain df.anno containing annotations"))
    }
    if (missing(colData)) {
        colData <- colData.Rds$df.anno
    } else {
        colData <- as.data.frame(colData)
        if (!("sample" %in% colnames(colData))) {
            colnames(colData)[1] <- "sample"
        }
        if (!all(colData$sample %in% colData.Rds$df.anno$sample)) {
            .log(paste("In makeSE():",
                "some samples in colData were not found in given path"),
                "message")
            colData <- colData[colData$sample %in% colData.Rds$df.anno$sample, ]
        }
    }
    return(colData)
}

# Converts charactor vectors to factors, removes columns with all NA's
.makeSE_colData_clean <- function(colData) {
    remove_na <- NULL
    if (ncol(colData) > 1) {
        for (i in seq(2, ncol(colData))) {
            if (is(colData[, i], "character")) {
                colData[, i] <- factor(unlist(colData[, i]))
            } else if (is(colData[, i], "logical")) {
                colData[, i] <- factor(unlist(
                    ifelse(colData[, i], "TRUE", "FALSE")))
            } else if (all(is.na(unlist(colData[, i])))) {
                remove_na <- append(remove_na, i)
            }
        }
    }
    if (!is.null(remove_na)) {
        colData <- colData[, -remove_na]
    }
    return(colData)
}

# Loads a NxtSE RDS
.makeSE_load_NxtSE <- function(filepath) {
    se <- readRDS(filepath)
    path <- dirname(filepath)
    se@assays <- .collateData_expand_assay_paths(se@assays, path)
    se@metadata[["Up_Inc"]] <- .collateData_expand_assay_path(
        se@metadata[["Up_Inc"]], path)
    se@metadata[["Down_Inc"]] <- .collateData_expand_assay_path(
        se@metadata[["Down_Inc"]], path)
    se@metadata[["Up_Exc"]] <- .collateData_expand_assay_path(
        se@metadata[["Up_Exc"]], path)
    se@metadata[["Down_Exc"]] <- .collateData_expand_assay_path(
        se@metadata[["Down_Exc"]], path)
    return(se)
}

# Iterates through IRFinder introns; removes overlapping minor introns
.makeSE_iterate_IR <- function(se, collate_path) {

    junc_PSI <- HDF5Array(file.path(normalizePath(collate_path),
        "data.h5"), "junc_PSI")[, colnames(se), drop = FALSE]

    se.IR <- se[rowData(se)$EventType == "IR", , drop = FALSE]
    se.coords <- rowData(se.IR)$EventRegion[
        rowData(se.IR)$EventRegion %in% rownames(junc_PSI)]
    se.coords.gr = coord2GR(se.coords)
    names(se.coords.gr) = se.coords
    
    if (length(se.coords.gr) > 0) {
        .log(paste("Iterating through IR events to determine introns",
            "of main isoforms"), type = "message")
        include <- .makeSE_iterate_IR_select_events(se.coords.gr, junc_PSI)
        se.coords.final <- se.coords.gr[include]
        se.coords.excluded <- se.coords.gr[!include]

        # Iteration to find events not overlapping with se.IR.final
        include <- .makeSE_iterate_IR_retrieve_excluded_introns(
            se.coords.final, se.coords.excluded)
        iteration <- 0
        while (length(include) > 0 & length(se.coords.final) > 0) {
            iteration <- iteration + 1
            .log(paste("Iteration", iteration), type = "message")
            dash_progress(paste("Iteration", iteration), 8)
            se.coords.excluded <- se.coords.excluded[include]

            include <- .makeSE_iterate_IR_select_events(
                se.coords.excluded, junc_PSI)

            if (length(include) > 0 && !all(include)) {
                se.coords.final <- c(se.coords.final,
                    se.coords.excluded[include])
                se.coords.excluded <-
                    se.coords.excluded[!include]
                include <- .makeSE_iterate_IR_retrieve_excluded_introns(
                    se.coords.final, se.coords.excluded)
            } else if (length(include) > 0) {
                se.coords.final <- c(se.coords.final,
                    se.coords.excluded)
                include <- c()
            } else {
                include <- c()
            }
        }

        se <- se[c(
            which(rowData(se.IR)$EventRegion %in% names(se.coords.final)),
            which(rowData(se)$EventType != "IR")
        ), ]
    }
    return(se)
}

# Selects introns of major isoforms
.makeSE_iterate_IR_select_events <- function(se.coords.gr, junc_PSI) {
    if(length(se.coords.gr) == 0) return(logical(0))
    if(length(se.coords.gr) == 1) return(TRUE)
    
    gr <- se.coords.gr
    gr.reduced <- reduce(gr)

    OL <- findOverlaps(gr, gr.reduced)
    junc_PSI.group <- as.data.table(junc_PSI[names(se.coords.gr), , drop = FALSE])
    junc_PSI.group$means <- rowMeans(junc_PSI.group)
    junc_PSI.group$group <- to(OL)
    junc_PSI.group[, c("max_means") := max(get("means")),
        by = "group"]
    return(junc_PSI.group$means == junc_PSI.group$max_means)
}

# Find excluded introns that does not overlap with given selection of introns
.makeSE_iterate_IR_retrieve_excluded_introns <- function(
        se.coords.final, se.coords.excluded) {
        
    if (length(se.coords.excluded) > 0) {
        if(length(se.coords.final) == 0) {
            return(rep(TRUE, length(se.coords.excluded)))
        }
        final.gr <- se.coords.final
        excluded.gr <- se.coords.excluded

        OL <- findOverlaps(excluded.gr, final.gr)
        include <- which(!(
            seq_len(length(excluded.gr))) %in% sort(unique(from(OL))))
    } else {
        include <- c()
    }
    return(include)
}
