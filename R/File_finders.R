#' Convenience Function to (recursively) find all files in a folder.
#'
#' Often, files e.g. alignment BAM files,
#' or [processBAM] output files, are stored in a single folder under some 
#' directory structure.
#' They can be grouped by being in common directory or having common names.
#' Often, their sample names can be gleaned by these common names or the names
#' of the folders in which they are contained.
#' This function (recursively) finds all files and
#' extracts sample names assuming either the files are named by sample names
#' (`level = 0`), or that their names can be derived from the
#' parent folder (`level = 1`). Higher `level` also work (e.g. `level = 2`)
#' mean the parent folder of the parent folder of the file is named by sample
#' names. See details section below.
#'
#' @details
#'
#' In BAM files, often the parent directory denotes their sample names. In this
#' case, use `level = 1` to automatically annotate the sample names using
#' `findBAMS()`.
#'
#' [processBAM] outputs two files per BAM processed. These are named by the 
#' given sample names. The text output is named "sample1.txt.gz", and the COV
#' file is named "sample1.cov", where `sample1` is the name of the sample. These
#' files can be organised / tabulated using the function `findSpliceWizOutput`.
#' The generic function `findSamples` will organise the [processBAM] text output
#' files but exclude the COV files. Use the latter as the `Experiment` in
#' [collateData] if one decides to collate an experiment without linked COV
#' files, for portability reasons.
#'
#' @param sample_path The path in which to recursively search for files
#'   that match the given `suffix`
#' @param suffix A vector of or or more strings that specifies the file suffix
#'   (e.g. '.bam' denotes BAM files, whereas ".txt.gz" denotes gzipped txt
#'   files).
#' @param level Whether sample names can be found in the file names themselves
#'   (level = 0), or their parent directory (level = 1). Potentially parent
#'   of parent directory (level = 2). Support max level <= 3 (for sanity).
#' @return A multi-column data frame with the first column containing
#'   the sample name, and subsequent columns being the file paths with suffix
#'   as determined by `suffix`.
#' @examples
#' # Retrieve all BAM files in a given folder, named by sample names
#' bam_path <- tempdir()
#' example_bams(path = bam_path)
#' df.bams <- findSamples(sample_path = bam_path,
#'   suffix = ".bam", level = 0)
#' # equivalent to:
#' df.bams <- findBAMS(bam_path, level = 0)
#'
#' # Retrieve all processBAM() output files in a given folder,
#' # named by sample names
#'
#' expr <- findSpliceWizOutput(file.path(tempdir(), "SpliceWiz_Output"))
#'
#' @name findSamples
#' @md
NULL

#' @describeIn findSamples Finds all files with the given suffix pattern.
#' Annotates sample names based on file or parent folder names.
#' @export
findSamples <- function(sample_path, suffix = ".txt.gz", level = 0) {
    if (length(suffix) == 0)
        .log(paste("In findSamples(),",
            "suffix must be of length greater than zero"))

    if (!dir.exists(sample_path))
        .log(paste("In findSamples(),",
            sample_path, "- given path does not exist"))

    if (length(level) > 1) .log("level must be a numeral of length 1")
    if (!is.numeric(level)) .log("In findSamples(), level must be numeric")
    if (level < 0) .log("In findSamples(), level must be non-negative")

    level <- floor(level)
    suffix_name <- "path"
    if (length(suffix) > 1) suffix_name <- suffix

    DT.list <- list()
    for (i in seq_len(length(suffix))) {
        pattern <- paste0("\\", suffix[i], "$")
        files_found <- list.files(pattern = pattern,
            path = normalizePath(sample_path),
            full.names = TRUE, recursive = TRUE)
        if (length(files_found) > 0) {
            DT <- data.table(sample = gsub(pattern, "",
                    files_found),
                path = files_found)
            lvl <- level
            while (lvl > 0) {
                DT$sample <- dirname(DT$sample)
                if (any(DT$sample %in% c(".", "/", "~")))
                    .log(paste("Sample name points to", DT$sample,
                        "- please check level of sample names"))
                lvl <- lvl - 1
            }
            DT$sample <- basename(DT$sample)
            colnames(DT)[2] <- suffix_name[i]
            setkeyv(DT, "sample")
            DT.list[[i]] <- DT
        } else {
            DT.list[[i]] <- NULL
        }
    }
    if (length(DT.list) <= 1) return(as.data.frame(DT.list))

    final <- DT.list[[1]]
    if (length(suffix) > 1) {
        # Check identity of sorted sample names
        samples <- DT.list[[1]]$sample
        for (i in seq(2, length(suffix))) {
            if (!is.null(DT.list[[i]]) &&
                    identical(samples, DT.list[[i]]$sample)) {
                final <- cbind(final, DT.list[[i]][, 2, with = FALSE])
            }
        }
    }
    cols <- c("sample", suffix_name[suffix_name %in% colnames(final)])
    return(as.data.frame(final[, cols, with = FALSE]))
}

#' @describeIn findSamples Use findSamples() to return all BAM files in a
#' given folder
#' @export
findBAMS <- function(sample_path, level = 0) {
    return(findSamples(sample_path, ".bam", level = level))
}

#' @describeIn findSamples Use findSamples() to return all processBAM output
#' files in a given folder, including COV files
#' @export
findSpliceWizOutput <- function(sample_path, level = 0) {
    DT <- findSamples(sample_path, c(".txt.gz", ".cov"), level = level)
    colnames(DT) <- c("sample", "sw_file", "cov_file")
    return(DT)
}
