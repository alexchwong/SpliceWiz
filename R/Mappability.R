#' Calculate low mappability genomic regions
#'
#' @description
#' These functions empirically calculate low-mappability (Mappability Exclusion)
#' regions using the given genome FASTA file. A splice-aware alignment software
#' capable of aligning reads to the genome is required.
#' See details and examples below.
#'
#' @details Creating a Mappability Exclusion BED file is a three-step process.
#' \cr
#' * First, using `Mappability_GenReads()`,
#' synthetic reads are systematically generated using the given genome contained
#' within `reference_path`.
#' * Second, an aligner
#' such as STAR (preferably the same aligner used for the subsequent RNA-seq
#' experiment) is required to align these reads to the source genome. Poorly
#' mapped regions of the genome will be reflected by regions of low coverage
#' depth.
#' * Finally, the BAM file containing the aligned reads is analysed
#' using `Mappability_CalculateExclusions()`, to identify
#' low-mappability regions to compile the Mappability Exclusion BED file.
#'
#' It is recommended to leave all parameters to their default settings. Regular
#' users should only specify `reference_path`, `aligned_bam` and `n_threads`,
#' as required.
#'
#' NB: [STAR_Mappability] runs all 3 steps required, using the `STAR` aligner.
#' This only works in systems where `STAR` is installed.
#'
#' NB2: In systems where `STAR` is not available, consider using HISAT2 or
#' Rsubread. A working example using Rsubread is shown below.
#'
#' Additionally, SpliceWiz retrieves pre-generated Mappability exclusion regions
#' suitable for use in generating references based on hg38,
#' hg19, mm10 and mm9 genomes. These were generated empirically. Synthetic 70-nt
#' reads, with start distances 10-nt apart, were systematically generated from
#' the genome. These reads were aligned to the same genome using the STAR 
#' aligner. Then, the BAM file read coverage was assessed.
#' Whereas mappable regions are expected to be covered with 7 reads,
#' low mappability regions are defined as regions covered with 4 or fewer
#' reads.
#'
#' @param genome_type Either one of `hg38`, `hg19`, `mm10` or `mm9`
#' @param as_type (Default "GRanges") Whether to return the Mappability 
#'   exclusion data as a GRanges object `"GRanges"`, or a BED `"bed"` or gzipped
#'   BED `"bed.gz"` copied locally to the given directory `path`.
#' @param path (Default = tempdir()) The desired destination path in which to 
#'   place a copy of the files. The directory does not need to exist but its 
#'   parent directory does.
#' @param overwrite (Default = `FALSE`)
#'   Whether or not to overwrite files if they already exist in the given path. 
#' @param offline (Default = `FALSE`)
#'   Whether or not to work in offline mode. This may be suitable
#'   if these functions have been previously run and the user wishes to run
#'   these functions without fetching online hub resources. Default = FALSE
#' @param reference_path The directory of the reference prepared by
#'   `GetReferenceResource()`
#' @param read_len The nucleotide length of the synthetic reads
#' @param read_stride The nucleotide distance between adjacent synthetic reads
#' @param error_pos The position of the procedurally-generated nucleotide error
#'   from the start of each synthetic reads
#' @param verbose Whether additional status messages are shown
#' @param alt_fasta_file (Optional) The path to the user-supplied genome fasta
#'   file, if different to that found inside the `resource` subdirectory of the
#'   `reference_path`. If `GetReferenceResource()` has already been run,
#'   this parameter should be omitted.
#' @param aligned_bam The BAM file of alignment of the synthetic reads generated
#'   by `Mappability_GenReads()`. Users should use a genome splice-aware
#'   aligner, preferably the same aligner used to align the samples in their
#'   experiment.
#' @param threshold Genomic regions with this alignment read depth (or below)
#'   in the aligned synthetic read BAM are defined as low
#'   mappability regions.
#' @param n_threads The number of threads used to calculate mappability
#'   exclusion regions from aligned bam file of synthetic reads.
#' @return
#' * For `Mappability_GenReads`: writes `Reads.fa` to the `Mappability`
#'   subdirectory inside the given `reference_path`.
#' * For `Mappability_CalculateExclusions`: writes a gzipped BED file
#'   named `MappabilityExclusion.bed.gz` to the `Mappability` subdirectory
#'   inside `reference_path`.
#'   This BED file is automatically used by `BuildReference()` if
#'   `MappabilityRef` is not specified.
#' @examples
#'
#' # (0) (Optional) Retrieve ready-to-use Mappability Exclusion reference
#' #                from AnnotationHub:
#'
#' # Returns the Mappability exclusion for hg38 directly as GRanges object
#'
#' hg38.MapExcl.gr <- get_mappability_exclusion(
#'     genome_type = "hg38", 
#'     as_type = "GRanges"
#' ) 
#' 
#' # returns the location of the Mappability exclusion gzipped BED for hg38
#'
#' gzippedBEDpath <- get_mappability_exclusion(
#'     genome_type = "hg38", 
#'     as_type = "bed.gz",
#'     path = tempdir()
#' ) 
#'
#' # (1a) Creates genome resource files
#'
#' ref_path <- file.path(tempdir(), "Reference")
#'
#' GetReferenceResource(
#'     reference_path = ref_path,
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#'
#' # (1b) Systematically generate reads based on the NxtIRF example genome:
#'
#' Mappability_GenReads(
#'     reference_path = ref_path
#' )
#' \dontrun{
#'
#' # (2) Align the generated reads using Rsubread:
#'
#' # (2a) Build the Rsubread genome index:
#'
#' setwd(ref_path)
#' Rsubread::buildindex(basename = "./reference_index",
#'     reference = chrZ_genome())
#'
#' # (2b) Align the synthetic reads using Rsubread::subjunc()
#'
#' Rsubread::subjunc(
#'     index = "./reference_index",
#'     readfile1 = file.path(ref_path, "Mappability", "Reads.fa"),
#'     output_file = file.path(ref_path, "Mappability", "AlignedReads.bam"),
#'     useAnnotation = TRUE,
#'     annot.ext = chrZ_gtf(),
#'     isGTF = TRUE
#' )
#'
#' # (3) Analyse the aligned reads in the BAM file for low-mappability regions:
#'
#' Mappability_CalculateExclusions(
#'     reference_path = ref_path,
#'     aligned_bam = file.path(ref_path, "Mappability", "AlignedReads.bam")
#' )
#'
#' # (4) Build the SpliceWiz reference using the calculated Mappability 
#' #     Exclusions
#'
#' BuildReference(ref_path)
#'
#' # NB the default is to search for the BED file generated by
#' # `Mappability_CalculateExclusions()` in the given reference_path
#' }
#' @name Mappability-methods
#' @aliases
#' Mappability_GenReads
#' Mappability_CalculateExclusions
#' @seealso [BuildReference]
#' @md
NULL

#' @describeIn Mappability-methods Fetches Mappability Exclusion reference
#' from ExperimentHub and places a copy in the given path; 
#' returns the location of this Mappability exclusion resource file
#' @export
get_mappability_exclusion <- function(
        genome_type = c("hg38", "hg19", "mm10", "mm9"),
        as_type = c("GRanges", "bed", "bed.gz"),
        path = tempdir(), overwrite = FALSE, offline = FALSE
) {
    genome_type <- match.arg(genome_type)
    if(genome_type == "") 
        stop("genome_type must be one of `hg38`, `hg19`, `mm10`, or `mm9`")
    as_type <- match.arg(as_type)
    if(as_type == "") 
        stop("as_type must be one of `GRanges`, `bed`, or `bed.gz`")
    if(as_type != "GRanges") .find_and_create_dir(path)

    title <- paste("NxtIRF", "mappability", genome_type, sep="/")
    destfile <- sprintf(file.path(path, "%s.MappabilityExclusion.bed"),
        genome_type)
    if(file.exists(paste0(destfile, ".gz")) & !overwrite & as_type != "bed.gz") 
        return(paste0(destfile, ".gz"))
    if(file.exists(destfile) & !overwrite & as_type != "bed") return(destfile)
    
    hubobj <- NULL
    # Check cache:
    if(!is.na(.query_local_cache(title))) {
        gr <- readRDS(.query_local_cache(title))
    } else {
        tryCatch({
            hubobj <- ExperimentHub(localHub = offline)
        }, error = function(e) {
            hubobj <- NULL
        })
        if(is.null(hubobj)) {
            message(
                "Failed establishing ExperimentHub connection. ",
                "Run ExperimentHub() to reproduce error msg"
            )
            return(NULL)
        }
        record_name <- names(hubobj[hubobj$title == title])
        if(length(record_name) < 1) {
            stopmsg <- paste("Mappability record not found -", genome_type,
                ifelse(offline, "- Perhaps try again in `online` mode.",
                paste("- Ensure ExperimentHub package is",
                    "updated to the latest version")))
            stop(stopmsg)
        } else if(length(record_name) > 1) {
            stopmsg <- paste("Multiple mappability records found -", 
                genome_type,
                "- please update SpliceWiz to latest version")
            stop(stopmsg)
        }
        tryCatch({
            cache_loc <- cache(hubobj[hubobj$title == title])
        }, error = function(e) {
            cache_loc <- ""
        })
        if(!file.exists(cache_loc)) {
            message("Downloading mappability from ExperimentHub failed")
            return(NULL)
        }
        .add_file_to_local_cache(cache_loc, title)
        gr <- hubobj[[record_name]]  # GRanges object from Rds
    }
    
    if(as_type == "GRanges") return(gr)
    if(!file.exists(destfile) | overwrite) {
        if(file.exists(destfile)) file.remove(destfile)
        rtracklayer::export(gr, destfile, "bed")
    }
    if(!file.exists(destfile)) {
        message("rtracklayer BED export failed for - ", genome_type)
        return(NULL)
    }
    if(as_type == "bed") return(destfile)
    if(file.exists(paste0(destfile, ".gz")))
        file.remove(paste0(destfile, ".gz"))
    R.utils::gzip(destfile)
    return(paste0(destfile, ".gz"))
}

#' @describeIn Mappability-methods Generates synthetic reads from a
#' genome FASTA file, for mappability calculations.
#' @export
Mappability_GenReads <- function(reference_path,
        read_len = 70, read_stride = 10, error_pos = 35,
        verbose = TRUE, alt_fasta_file) {
    .gmr_check_params(read_len, read_stride, error_pos + 1)
    if (missing(alt_fasta_file)) {
        alt_fasta_file <- .STAR_get_FASTA(reference_path)
        if (!file.exists(alt_fasta_file)) .log(paste("In Mappability_GenReads,",
            "failed to generate genome fasta file from given reference"))
    } else if (!file.exists(alt_fasta_file)) {
        .log(paste("In Mappability_GenReads,",
            "given fasta file", alt_fasta_file, "not found"))
    }
    .validate_path(file.path(normalizePath(reference_path), "Mappability"))
    # Run map read generator:
    outfile <- file.path(normalizePath(reference_path),
        "Mappability", "Reads.fa")
    .log(paste("Generating synthetic reads, saving to", outfile), "message")
    .run_IRFinder_GenerateMapReads(
        normalizePath(alt_fasta_file), outfile,
        read_len, read_stride, error_pos + 1
    )
    .STAR_clean_temp_FASTA_GTF(reference_path)
}

#' @describeIn Mappability-methods Generate a BED file defining
#' low mappability regions, using reads generated by
#' \code{Mappability_GenReads()}, aligned to the genome.
#' @export
Mappability_CalculateExclusions <- function(reference_path,
        aligned_bam = file.path(reference_path, "Mappability",
            "Aligned.out.bam"),
        threshold = 4, n_threads = 1) {
    if (!file.exists(aligned_bam))
        .log(paste("In Mappability_CalculateExclusions(),",
            aligned_bam, "BAM file does not exist"))

    .validate_path(file.path(normalizePath(reference_path), "Mappability"))
    output_file <- file.path(normalizePath(reference_path), "Mappability",
        "MappabilityExclusion.bed")

    .log(paste("Calculating Mappability Exclusion regions from:",
        aligned_bam), type = "message")
    .run_IRFinder_MapExclusionRegions(
        bamfile = normalizePath(aligned_bam),
        output_file = output_file,
        threshold = threshold,
        n_threads = n_threads
    )
}

# Checks whether given parameters are valid
.gmr_check_params <- function(read_len, read_stride, error_pos) {
    if (!is.numeric(read_len) || read_len < 30) {
        .log(paste("In Mappability_GenReads,",
            "read_len must be numerical and at least 30"))
    }
    if (!is.numeric(read_stride) || read_stride > read_len) {
        .log(paste("In Mappability_GenReads,",
            "read_stride must be numerical and less than read_len"))
    }
    if (!is.numeric(error_pos) || error_pos > read_len) {
        .log(paste("In Mappability_GenReads,",
            "error_pos must be numerical and less than read_len"))
    }
}

# Wrappers to native Rcpp functions:

.run_IRFinder_GenerateMapReads <- function(genome.fa = "", out.fa,
    read_len = 70, read_stride = 10, error_pos = 36) {
    return(
        c_GenerateMappabilityReads(normalizePath(genome.fa),
            file.path(normalizePath(dirname(out.fa)), basename(out.fa)),
            read_len = read_len,
            read_stride = read_stride,
            error_pos = error_pos)
    )
}

.run_IRFinder_MapExclusionRegions <- function(bamfile = "", output_file,
        threshold = 4, includeCov = FALSE, n_threads = 1) {
    s_bam <- normalizePath(bamfile)

    c_GenerateMappabilityRegions(s_bam,
        output_file,
        threshold = threshold,
        includeCov = includeCov,
        verbose = TRUE, n_threads = n_threads
    )
    # check file is actually made; then gzip it
    if (file.exists(paste0(output_file, ".txt"))) {
        R.utils::gzip(filename = paste0(output_file, ".txt"),
            destname = paste0(output_file, ".gz"))
    } else {
        .log(paste(paste0(output_file, ".txt"), "was not produced"))
    }
}
