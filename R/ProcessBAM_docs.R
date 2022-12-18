#' Runs the OpenMP/C++ based SpliceWiz algorithm
#'
#' These function calls the SpliceWiz C++ routine on one or more BAM files.
#' \cr\cr The routine is an improved version over the original IRFinder, with
#' OpenMP-based multi-threading and the production of compact "COV" files to
#' record alignment coverage. A SpliceWiz reference built using 
#' [Build-Reference-methods] is required.\cr\cr
#' After `processBAM()` is run, users should call
#' [collateData] to collate individual outputs into an experiment / dataset.
#' \cr\cr
#' BAM2COV creates COV files from BAM files without running `processBAM()`.
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
#' @param verbose (default `FALSE`) Set to `TRUE` to allow `processBAM()` to 
#'   output progress bars and messages
#' @param multiRead (default `FALSE`) Whether SpliceWiz/ompBAM should use one
#'   (set to `FALSE`) or all available threads (set to `TRUE`) to read BAM
#'   files from the storage drive. In SSD drives or high performance computing
#'   clusters, setting to `TRUE` may slightly improve performance, whereas if
#'   reading from disk is the speed bottleneck, the default setting `FALSE`
#'   should result in higher performance.
#' @return
#' Output will be saved to `output_path`. Output files 
#'   will be named using the given `sample_names`.
#' For `processBAM()`:
#'
#'   * sample.txt.gz: The main output file containing the quantitation
#'   of IR and splice junctions, as well as QC information\cr\cr
#'   * sample.cov: Contains coverage information in compressed binary. See
#'     [getCoverage]
#'   * main.FC.Rds: A single file containing gene counts for the whole dataset
#'   (only if `run_featureCounts == TRUE`)
#'
#' For `BAM2COV()`:
#'
#'   * sample.cov: Contains coverage information in compressed binary. See
#'     [getCoverage]
#'
#' @examples
#'
#' # Run BAM2COV, which only produces COV files but does not run `processBAM()`:
#'
#' bams <- SpliceWiz_example_bams()
#'
#' BAM2COV(bams$path, bams$sample,
#'   output_path = file.path(tempdir(), "SpliceWiz_Output"),
#'   n_threads = 2, overwrite = TRUE
#' )
#'
#' # Run processBAM(), which produces:
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
#' processBAM(bams$path, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "SpliceWiz_Output"),
#'   n_threads = 2
#' )
#' @seealso [Build-Reference-methods] [collateData] [isCOV]
#' @name processBAM
#' @md
NULL