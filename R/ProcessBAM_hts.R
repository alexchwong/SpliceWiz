#' Runs the OpenMP/C++ based SpliceWiz algorithm (htslib version)
#'
#' These function calls the SpliceWiz C++ routine on one or more BAM files.
#' \cr\cr The routine is an improved version over the original IRFinder, with
#' OpenMP-based multi-threading and the production of compact "COV" files to
#' record alignment coverage. A SpliceWiz reference built using 
#' [Build-Reference-methods] is required.\cr\cr
#' After `processBAM_hts()` is run, users should call
#' [collateData] to collate individual outputs into an experiment / dataset.
#' \cr\cr
#' BAM2COV creates COV files from BAM files without running `processBAM_hts()`.
#' \cr\cr See details for performance info.
#'
#' @details
#' Typical run-times for a 100-million paired-end alignment BAM file takes 10
#' minutes using a single core. Using 8 threads, the runtime is approximately
#' 2-5 minutes, depending on your system's file input / output speeds. 
#' Approximately 10 Gb of RAM is used when OpenMP is used. If OpenMP
#' is not used (see below), this memory usage is multiplied across the number
#' of processor threads (i.e. 40 Gb if `n_threads = 4`).
#'
#' OpenMP is natively available to Linux / Windows compilers, and OpenMP will
#' be used if `useOpenMP` is set to `TRUE`, using multiple threads to process
#' each BAM file. On Macs, if OpenMP is not available at compilation,
#' BiocParallel will be used, processing BAM files simultaneously,
#' with one BAM file per thread.
#'
#' @param bamfiles A vector containing file paths of 1 or more BAM files
#' @param sample_names The sample names of the given BAM files. Must
#'   be a vector of the same length as `bamfiles`
#' @param reference_path The directory containing the SpliceWiz reference
#' @param output_path The output directory of this function
#' @param n_threads (default `1`) The number of threads to use. See details.
#' @param useOpenMP (default `TRUE`) Whether to use OpenMP.
#'   If set to `FALSE`, BiocParallel will be used if `n_threads` is set
#' @param overwrite (default `FALSE`) If output files already exist,
#'   will not attempt to re-run. If `run_featureCounts` is `TRUE`, will not
#'   overwrite gene counts of previous run unless `overwrite` is `TRUE`.
#' @param run_featureCounts (default `FALSE`) Whether this function will run
#'   [Rsubread::featureCounts] on the BAM files after counting spliced reads.
#'   If so, the output will be
#'   saved to `"main.FC.Rds` in the `output_path` directory as a list object.
#' @param verbose (default `FALSE`) Set to `TRUE` to allow `processBAM_hts()` to 
#'   output progress bars and messages
#' @param multiRead (default `FALSE`) Whether to use multiple threads to read
#'   files. Set to `TRUE` to use multiple threads (improves performance if
#'   files are already cached into memory.
#' @param read_pool_size How many reads to decompress (htslib) prior to 
#'   processing of reads. Default `1000000`
#' @return `processBAM_hts()` output will be saved to `output_path`. Output files 
#'   will be named using the given sample names.
#'   * sample.txt.gz: The main output file containing the quantitation
#'   of IR and splice junctions, as well as QC information\cr\cr
#'   * sample.cov: Contains coverage information in compressed binary. See
#'     [getCoverage]
#'   * main.FC.Rds: A single file containing gene counts for the whole dataset
#'   (only if `run_featureCounts == TRUE`)
#' @examples
#'
#' # Run BAM2COV, which only produces COV files but does not run `processBAM_hts()`:
#'
#' bams <- SpliceWiz_example_bams()
#'
#' BAM2COV_hts(bams$path, bams$sample,
#'   output_path = file.path(tempdir(), "SpliceWiz_Output"),
#'   n_threads = 2, overwrite = TRUE
#' )
#'
#' # Run processBAM_hts(), which produces:
#' # - text output of intron coverage and spliced read counts
#' # - COV files which record read coverages
#'
#' example_ref <- file.path(tempdir(), "Reference")
#'
#' buildRef(
#'     reference_path = example_ref,
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#'
#' bams <- SpliceWiz_example_bams()
#'
#' processBAM_hts(bams$path, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "SpliceWiz_Output"),
#'   n_threads = 2
#' )
#' @seealso [Build-Reference-methods] [collateData] [isCOV]
#' @name processBAM_hts
#' @md
NULL

#' @describeIn processBAM_hts Converts BAM files to COV files without running
#'   `processBAM_hts()`
#' @export
BAM2COV_hts <- function(
        bamfiles = "./Unsorted.bam",
        sample_names = "sample1",
        output_path = "./cov_folder",
        n_threads = 1, useOpenMP = TRUE,
        overwrite = FALSE,
        verbose = FALSE,
        read_pool_size = 1000000
) {
    # Check args
    if (length(bamfiles) != length(sample_names)) 
        .log(paste("In BAM2COV_hts(),",
            "Number of BAM files and sample names must be the same"))
    if (length(sample_names) != length(unique(sample_names)))
        .log(paste("In BAM2COV_hts(), some sample names are not unique"))

    if (length(bamfiles) == 0) .log("bamfiles argument must not be empty")
    if (!all(file.exists(bamfiles)))
        .log(paste("In BAM2COV_hts(),", "some BAMs in bamfiles do not exist"))

    if (!dir.exists(dirname(output_path))) .log(paste("In BAM2COV_hts(),",
        dirname(output_path), " - path does not exist"))

    if (!dir.exists(output_path)) dir.create(output_path)

    # Check which output already exists; prevent overwrite
    s_output <- file.path(normalizePath(output_path), sample_names)
    if (!overwrite) {
        already_exist <- (
            file.exists(paste0(s_output, ".cov"))
        )
    } else {
        already_exist <- rep(FALSE, length(bamfiles))
    }

    # Call wrapper
    if (!all(already_exist)) {
        .run_BAM2COV_hts(
            bamfiles = bamfiles[!already_exist],
            output_file_prefixes = s_output[!already_exist],
            max_threads = n_threads, useOpenMP = useOpenMP,
            overwrite = overwrite,
            verbose = verbose, read_pool_size = read_pool_size
        )
    } else {
        .log("BAM2COV has already been run on given BAM files", "message")
    }

    s_output <- file.path(normalizePath(output_path), sample_names)
    if (!all(file.exists(paste0(s_output, ".cov"))))
        .log(paste("Some BAM2COV outputs could not be found.",
            "BAM2COV must have crashed"))
}

#' @describeIn processBAM_hts Processes BAM files. Requires a
#' SpliceWiz reference generated by buildRef()
#' @export
processBAM_hts <- function(
        bamfiles = "./Unsorted.bam",
        sample_names = "sample1",
        reference_path = "./Reference",
        output_path = "./SpliceWiz_Output",
        n_threads = 1, useOpenMP = TRUE,
        overwrite = FALSE,
        run_featureCounts = FALSE,
        verbose = FALSE,
        read_pool_size = 1000000
) {
    # Check args
    if (length(bamfiles) != length(sample_names)) .log(paste("In processBAM_hts(),",
        "Number of BAM files and sample names must be the same"))
    if (length(sample_names) != length(unique(sample_names)))
        .log(paste("In processBAM_hts(), some sample names are not unique"))

    if (length(bamfiles) == 0) .log("bamfiles argument must not be empty")
    if (!all(file.exists(bamfiles))) .log(paste("In processBAM_hts(),",
        "some BAMs in bamfiles do not exist"))

    if (!dir.exists(dirname(output_path))) .log(paste("In processBAM_hts(),",
        dirname(output_path), " - path does not exist"))

    if (!dir.exists(output_path)) dir.create(output_path)

    # Check which output already exists; prevent overwrite
    s_output <- file.path(normalizePath(output_path), sample_names)
    if (!overwrite) {
        already_exist <- (
            file.exists(paste0(s_output, ".txt.gz")) &
            file.exists(paste0(s_output, ".cov"))
        )
    } else {
        already_exist <- rep(FALSE, length(bamfiles))
    }

    # Call wrapper
    if (!all(already_exist)) {
        .run_processBAM_hts(
            reference_path = reference_path,
            bamfiles = bamfiles[!already_exist],
            output_files = s_output[!already_exist],
            max_threads = n_threads, useOpenMP = useOpenMP,
            overwrite_SpliceWiz_Output = overwrite,
            verbose = verbose, read_pool_size = read_pool_size
        )
    } else {
        .log("processBAM has already been run on given BAM files", "message")
    }

    s_output <- file.path(normalizePath(output_path), sample_names)
    if (!all(file.exists(paste0(s_output, ".txt.gz"))))
        .log(paste("Some processBAM outputs could not be found.",
            "processBAM must have crashed"))

    # Run featureCounts
    if (run_featureCounts) {
        .processBAM_run_featureCounts(
            reference_path, s_output,
            bamfiles, sample_names, n_threads, overwrite
        )
    }
}

# Wrapper to R/C++. Handles whether OpenMP or BiocParallel is used
.run_processBAM_hts <- function(
        reference_path = "./Reference",
        bamfiles = "Unsorted.bam",
        output_files = "./Sample",
        max_threads = max(parallel::detectCores(), 1),
        useOpenMP = TRUE,
        overwrite_SpliceWiz_Output = FALSE,
        verbose = TRUE, read_pool_size = 1000000
    ) {
    .validate_reference(reference_path) # Check valid SpliceWiz reference
    s_bam <- normalizePath(bamfiles) # Clean path name for C++
    s_ref <- normalizePath(reference_path) # Clean path name for C++

    # Check args
    .processBAM_validate_args(s_bam, max_threads, output_files)
    ref_file <- file.path(s_ref, "SpliceWiz.ref.gz")

    .log("Running SpliceWiz processBAM", "message")
    n_threads <- floor(max_threads)
    if (Has_OpenMP() > 0 & useOpenMP) {
        for (i in seq_len(length(s_bam))) {
            ret <- SpliceWizMain_hts(s_bam[i], ref_file, 
                output_files[i], verbose, n_threads, read_pool_size)
        }

    } else {
        # Use BiocParallel
        n_rounds <- ceiling(length(s_bam) / floor(max_threads))
        n_threads <- ceiling(length(s_bam) / n_rounds)

        BPPARAM_mod <- .validate_threads(n_threads, as_BPPARAM = TRUE)

        row_starts <- seq(1, by = n_threads, length.out = n_rounds)
        for (i in seq_len(n_rounds)) {
            selected_rows_subset <- seq(row_starts[i],
                min(length(s_bam), row_starts[i] + n_threads - 1)
            )
            BiocParallel::bplapply(selected_rows_subset,
                function(i, s_bam, reference_file,
                        output_files, verbose, overwrite) {
                    .processBAM_run_single_hts(s_bam[i], reference_file,
                        output_files[i], verbose, overwrite, read_pool_size)
                },
                s_bam = s_bam,
                reference_file = ref_file,
                output_files = output_files,
                verbose = verbose,
                overwrite = overwrite_SpliceWiz_Output,
                BPPARAM = BPPARAM_mod
            )
        }
    }

}

# BAM2COV wrapper to R/C++. Handles whether OpenMP or BiocParallel is used
.run_BAM2COV_hts <- function(
        bamfiles = "sample.bam",
        output_file_prefixes = "sample",
        max_threads = max(parallel::detectCores(), 1),
        useOpenMP = TRUE,
        overwrite = FALSE,
        verbose = TRUE,
        read_pool_size = 1000000
    ) {
    s_bam <- normalizePath(bamfiles) # Clean path name for C++
    # Check args
    .processBAM_validate_args(s_bam, max_threads, output_file_prefixes)

    .log("Running BAM2COV", "message")
    n_threads <- floor(max_threads)
    if (Has_OpenMP() > 0 & useOpenMP) {
        # Simple FOR loop:
        for (i in seq_len(length(s_bam))) {
            .BAM2COV_run_single_hts(s_bam[i], output_file_prefixes[i],
                n_threads, verbose, overwrite, read_pool_size)
        }
    } else {
        # Use BiocParallel
        n_rounds <- ceiling(length(s_bam) / floor(max_threads))
        n_threads <- ceiling(length(s_bam) / n_rounds)

        BPPARAM_mod <- .validate_threads(n_threads, as_BPPARAM = TRUE)

        row_starts <- seq(1, by = n_threads, length.out = n_rounds)
        for (i in seq_len(n_rounds)) {
            selected_rows_subset <- seq(row_starts[i],
                min(length(s_bam), row_starts[i] + n_threads - 1)
            )
            BiocParallel::bplapply(selected_rows_subset,
                function(i, s_bam, output_files, verbose, overwrite) {
                    .BAM2COV_run_single_hts(s_bam[i], output_files[i], 1,
                        verbose, overwrite, read_pool_size)
                },
                s_bam = s_bam,
                output_files = output_file_prefixes,
                verbose = verbose,
                overwrite = overwrite,
                BPPARAM = BPPARAM_mod
            )
        }
    }
}

# Call C++ on a single sample. Used for BiocParallel
.processBAM_run_single_hts <- function(
    bam, ref, out, verbose, overwrite, read_pool_size
) {
    file_gz <- paste0(out, ".txt.gz")
    file_cov <- paste0(out, ".cov")
    bam_short <- file.path(basename(dirname(bam)), basename(bam))
    if (overwrite ||
        !(file.exists(file_gz) | file.exists(file_cov))) {
        ret <- SpliceWizMain_hts(bam, ref, out, verbose, 1, read_pool_size)
        # Check SpliceWiz returns all files successfully
        if (ret != 0) {
            .log(paste(
                "SpliceWiz processBAM exited with errors"))
        } else if (!file.exists(file_gz)) {
            .log(paste(
                "SpliceWiz processBAM failed to produce", file_gz))
        } else if (!file.exists(file_cov)) {
            .log(paste(
                "SpliceWiz processBAM failed to produce", file_cov))
        } else {
            .log(paste("SpliceWiz processBAM processed", bam_short), "message")
        }
    } else {
        .log(paste("SpliceWiz processBAM output for", bam_short,
            "already exists, skipping..."), "message")
    }
}

# Call C++/BAM2COV on a single sample. Used for BiocParallel
.BAM2COV_run_single_hts <- function(
    bam, out, n_threads, verbose, overwrite, read_pool_size = 1000000
) {
    file_cov <- paste0(out, ".cov")
    bam_short <- file.path(basename(dirname(bam)), basename(bam))
    if (overwrite || !(file.exists(file_cov))) {
        ret <- c_BAM2COV_hts(bam, file_cov, verbose, n_threads, read_pool_size)
        # Check BAM2COV returns all files successfully
        if (ret != 0) {
            .log(paste(
                "BAM2COV exited with errors"))
        } else if (!file.exists(file_cov)) {
            .log(paste(
                "BAM2COV failed to produce", file_cov))
        } else {
            .log(paste("BAM2COV processed", bam_short), "message")
        }
    } else {
        .log(paste(file_cov,
            "already exists, skipping..."), "message")
    }
}

# Runs featureCounts on given BAM files, intended to be run after processBAM
# as processBAM determines the strandedness and paired-ness of the experiment
.processBAM_run_featureCounts <- function(
        reference_path, output_files,
        s_bam, s_names, n_threads, overwrite
) {
    .check_package_installed("Rsubread", "2.4.0")
    gtf_file <- Get_GTF_file(reference_path)

    # determine paired-ness, strandedness, assume all BAMS are the same
    data.list <- get_multi_DT_from_gz(
        normalizePath(paste0(output_files[1], ".txt.gz")),
        c("BAM", "Directionality")
    )
    stats <- data.list$BAM
    direct <- data.list$Directionality

    paired <- (stats$Value[3] == 0 & stats$Value[4] > 0) ||
        (stats$Value[3] > 0 && stats$Value[4] / stats$Value[3] / 1000)
    strand <- direct$Value[9]
    if (strand == -1) strand <- 2

    # Check which have already been run, do not run if overwrite = FALSE
    outfile <- file.path(dirname(output_files[1]), "main.FC.Rds")
    if (file.exists(outfile) & !overwrite) {
        res.old <- readRDS(outfile)
        need_to_do <- (!(s_names %in% res.old$targets))
    } else {
        need_to_do <- rep(TRUE, length(s_bam))
    }

    if (any(need_to_do)) {
        # Run FeatureCounts in bulk
        res <- Rsubread::featureCounts(
            s_bam[need_to_do],
            annot.ext = gtf_file,
            isGTFAnnotationFile = TRUE,
            strandSpecific = strand,
            isPairedEnd = paired,
            requireBothEndsMapped = paired,
            nthreads = n_threads
        )
        res$targets <- s_names[need_to_do]
        colnames(res$counts) <- s_names[need_to_do]
        colnames(res$stat)[-1] <- s_names[need_to_do]
        columns <- c("counts", "annotation", "targets", "stat")
        # Append to existing main.FC.Rds if exists, overwriting where necessary:
        if (file.exists(outfile)) {
            res.old <- readRDS(outfile)
            if (!all(columns %in% names(res))) {
                .log(paste(outfile,
                    "found but was not a valid SpliceWiz featureCounts",
                    "output; overwriting previous output"
                ), "warning")
            } else if (
                identical(res.old$annotation, res$annotation) &
                identical(res.old$stat$Status, res$stat$Status)
            ) {
                new_samples <- res$targets[!(res$targets %in% res.old$targets)]
                res$targets <- c(res.old$targets, new_samples)
                res$stat <- cbind(res.old$stat, res$stat[, new_samples])
                res$counts <- cbind(res.old$counts, res$counts[, new_samples])
            } else {
                .log(paste(
                    "featureCounts output not compatible with previous",
                    "output in", outfile, "; overwriting previous output"
                ), "warning")
            }
        }
        if (all(columns %in% names(res))) {
            saveRDS(res, outfile)
        } else {
            .log("Error encountered when running featureCounts")
        }
        .log(paste("featureCounts ran succesfully; saved to",
            outfile), "message")
    } else {
        .log("featureCounts has already been run on given BAM files", "message")
    }

}



# Validate arguments; return error if invalid
.processBAM_validate_args <- function(s_bam, max_threads, output_files) {
    if (!is.numeric(max_threads)) max_threads <- 1
    if (max_threads < 1) max_threads <- 1
    max_threads <- floor(max_threads)

    if (max_threads > 1 && max_threads > parallel::detectCores()) {
        .log(paste(
            max_threads, " threads is not allowed for this system"))
    }

    if (!all(file.exists(s_bam))) {
        .log(paste(
            paste(unique(s_bam[!file.exists(s_bam)]), collapse = ""),
            " - these BAM files were not found"))
    }

    if (!all(dir.exists(dirname(output_files)))) {
        .log(paste(
            paste(unique(dirname(
                    output_files[!dir.exists(dirname(output_files))])),
                collapse = ""),
            " - directories not found"))
    }

    if (!(length(s_bam) == length(output_files))) {
        .log("Number of output files and bam files must be the same")
    }
    return(TRUE)
}
