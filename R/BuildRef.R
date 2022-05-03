#' Builds reference files used by SpliceWiz
#'
#' @description
#' These function builds the reference required by the SpliceWiz engine, as well
#' as alternative splicing annotation data for SpliceWiz. See examples
#' below for guides to making the SpliceWiz reference.
#' @details
#' `getResources()` processes the files, downloads resources from
#' web links or from `AnnotationHub()`, and saves a local copy in the "resource"
#' subdirectory within the given `reference_path`. Resources are retrieved via
#' either:
#' 1. User-supplied FASTA and GTF file. This can be a file path, or a web link
#'   (e.g. 'http://', 'https://' or 'ftp://'). Use `fasta` and `gtf`
#'    to specify the files or web paths to use.
#' 2. AnnotationHub genome and gene annotation (Ensembl): supply the names of
#'    the genome sequence and gene annotations to `fasta` and `gtf`.
#'
#' `buildRef()` will first run `getResources()` if resources are
#' not yet saved locally (i.e. `getResources()` is not already run).
#' Then, it creates the SpliceWiz references. Typical run-times are
#' 5 to 10 minutes for human and mouse genomes (after resources are downloaded).
#'
#' NB: the parameters `fasta` and `gtf` can be omitted in `buildRef()` if
#' `getResources()` is already run.
#'
#' Typical usage involves running `buildRef()` for human and mouse genomes
#' and specifying the `genome_type` to use the default `MappabilityRef` and
#' `nonPolyARef` files for the specified genome. For non-human non-mouse
#' genomes, use one of the following alternatives:
#' * Create the SpliceWiz reference without using Mappability Exclusion regions. 
#'     To do this, simply run `buildRef()` and omit `MappabilityRef`. This is
#'     acceptable assuming the introns assessed are short and do not contain
#'     intronic repeats
#' * Calculating Mappability Exclusion regions using the STAR aligner,
#'     and building the SpliceWiz reference. This can be done using the
#'     `buildFullRef()` function, on systems where `STAR` is installed
#' * Instead of using the STAR aligner, any genome splice-aware aligner could be
#'     used. See [Mappability-methods] for details. After producing the
#'     `MappabilityExclusion.bed.gz` file (in the `Mappability` subfolder), run
#'     `buildRef()` using this file (or simply leave it blank).
#'
#' BED files are tab-separated text files containing 3 unnamed columns
#' specifying chromosome, start and end coordinates. To view an example BED
#' file, open the file specified in the path returned by
#' `getNonPolyARef("hg38")`
#'
#' See examples below for common use cases.
#'
#' @param reference_path (REQUIRED) The directory path to store the generated
#'   reference files
#' @param fasta The file path or web link to the user-supplied genome
#'   FASTA file. Alternatively, the name of the AnnotationHub record containing
#'   the genome resource. May be omitted if `getResources()` has already
#'   been run using the same `reference_path`.
#' @param gtf The file path or web link  to the user-supplied transcript
#'   GTF file (or gzipped GTF file). Alternatively, the name of the
#'   AnnotationHub record containing the transcript GTF file. May be omitted if
#'   `getResources()` has already been run using the same
#'   `reference_path`.
#' @param overwrite (default `FALSE`) For `getResources()`: if the
#'   genome FASTA and gene annotation GTF files already exist in the `resource`
#'   subdirectory, it will not be overwritten. For `buildRef()` and
#'   `buildFullRef()`: the SpliceWiz reference will not be overwritten
#'   if one already exist. A reference is considered to exist if
#'   the file `SpliceWiz.ref.gz` is present inside `reference_path`.
#' @param force_download (default `FALSE`) When online resources are retrieved,
#'   a local copy is stored in the `SpliceWiz` BiocFileCache. Subsequent calls
#'   to the web resource will fetch the local copy. Set `force_download` to
#'   `TRUE` will force the resource to be downloaded from the web. Set this to
#'   `TRUE` only if the web resource has been updated since the last retrieval.
#' @param chromosome_aliases (Highly optional) A 2-column data frame containing
#'   chromosome name conversions. If this is set, allows [processBAM] to parse
#'   BAM alignments to a genome whose chromosomes are named
#'   differently to the reference genome. The most common scenario is where
#'   Ensembl genome typically use chromosomes "1", "2", ..., "X", "Y", whereas
#'   UCSC/Gencode genome use "chr1", "chr2", ..., "chrX", "chrY". See example
#'   below. Refer to <https://github.com/dpryan79/ChromosomeMappings> for a
#'   list of chromosome alias resources.
#' @param genome_type Allows `buildRef()` to select default
#'   `nonPolyARef` and `MappabilityRef` for selected genomes. Allowed options
#'   are: `hg38`, `hg19`, `mm10`, and `mm9`.
#' @param nonPolyARef (Optional) A BED file of regions defining known
#'   non-polyadenylated transcripts. This file is used for QC analysis 
#'   to measure Poly-A enrichment quality of samples.
#'   If omitted, and `genome_type` is defined, the default for the specified
#'   genome will be used.
#' @param MappabilityRef (Optional) A BED file of low mappability regions due to
#'   repeat elements in the genome. If omitted, the file generated by
#'   [calculateMappability()] will be used where available, and if
#'   this is not, the default file for the specified `genome_type` will be used.
#'   If `genome_type` is not specified, `MappabilityRef` is not used.
#'   See details.
#' @param BlacklistRef A BED file of regions to be otherwise excluded from IR
#'   analysis. If omitted, a blacklist is not used (this is the default).
#' @param useExtendedTranscripts (default `TRUE`) Should non-protein-coding
#'   transcripts such as anti-sense and lincRNA transcripts be included in
#'   searching for IR / AS events? Setting `FALSE` (vanilla IRFinder) will
#'   exclude transcripts other than `protein_coding` and
#'   `processed_transcript` transcripts from IR analysis.
#' @param lowMemoryMode (default `TRUE`) By default, SpliceWiz converts FASTA
#'   files to TwoBit, then uses the TwoBit file to fetch genome sequences. In
#'   most cases, this method uses less memory and is faster, but can be very
#'   slow on some systems. Set this option to `FALSE` (which will convert the
#'   TwoBit file back to FASTA) if you experience
#'   very slow genome fetching (e.g. when annotating splice motifs).
#' @param n_threads The number of threads used to generate the STAR reference
#'   and mappability calculations. Multi-threading is not used for SpliceWiz
#'   reference generation (but multiple cores are utilised in data-table
#'   and fst file processing automatically, where available). See [STAR-methods]
#' @param use_STAR_mappability (default FALSE) In `buildFullRef()`,
#'   whether to run [STAR_mappability] to calculate low-mappability regions.
#'   We recommend setting this to `FALSE` for the common genomes
#'   (human and mouse), and to `TRUE` for genomes not supported by
#'   `genome_type`. When set to false, the MappabilityExclusion default file
#'   corresponding to `genome_type` will automatically be used.
#' @return
#' For `getResources`: creates the following local resources:
#' * `reference_path/resource/genome.2bit`: Local copy of the genome sequences
#'   as a TwoBitFile.
#' * `reference_path/resource/transcripts.gtf.gz`: Local copy of the gene
#'   annotation as a gzip-compressed file.
#' For `buildRef` and `buildFullRef`: creates a SpliceWiz reference
#'   which is written to the given directory specified by `reference_path`.
#'   Files created includes:
#' * `reference_path/settings.Rds`: An RDS file containing parameters used
#'   to generate the SpliceWiz reference
#' * `reference_path/SpliceWiz.ref.gz`: A gzipped text file containing collated
#'   SpliceWiz reference files. This file is used by [processBAM]
#' * `reference_path/fst/`: Contains fst files for subsequent easy access to
#'   SpliceWiz generated references
#' * `reference_path/cov_data.Rds`: An RDS file containing data required to
#'    visualise genome / transcript tracks.
#'
#' `buildFullRef` also creates a `STAR` reference located in the `STAR`
#'   subdirectory inside the designated `reference_path`
#'
#' For `getNonPolyARef`: Returns the file path to the BED file for
#' the nonPolyA loci for the specified genome.
#'
#' @examples
#' # Quick runnable example: generate a reference using SpliceWiz's example genome
#'
#' example_ref <- file.path(tempdir(), "Reference")
#' getResources(
#'     reference_path = example_ref,
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#' buildRef(
#'     reference_path = example_ref
#' )
#'
#' # NB: the above is equivalent to:
#'
#' example_ref <- file.path(tempdir(), "Reference")
#' buildRef(
#'     reference_path = example_ref,
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#'
#' # Get the path to the Non-PolyA BED file for hg19
#'
#' getNonPolyARef("hg19")
#' \dontrun{
#'
#' ### Long examples ###
#'
#' # Generate a SpliceWiz reference from user supplied FASTA and GTF files for a
#' # hg38-based genome:
#'
#' buildRef(
#'     reference_path = "./Reference_user",
#'     fasta = "genome.fa", gtf = "transcripts.gtf",
#'     genome_type = "hg38"
#' )
#'
#' # NB: Setting `genome_type = hg38`, will automatically use default
#' # nonPolyARef and MappabilityRef for `hg38`
#'
#' # Reference generation from Ensembl's FTP links:
#'
#' FTP <- "ftp://ftp.ensembl.org/pub/release-94/"
#' buildRef(
#'     reference_path = "./Reference_FTP",
#'     fasta = paste0(FTP, "fasta/homo_sapiens/dna/",
#'         "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
#'     gtf = paste0(FTP, "gtf/homo_sapiens/",
#'         "Homo_sapiens.GRCh38.94.chr.gtf.gz"),
#'     genome_type = "hg38"
#' )
#'
#' # Get AnnotationHub record names for Ensembl release-94:
#'
#' # First, search for the relevant AnnotationHub record names:
#'
#' ah <- AnnotationHub::AnnotationHub()
#' AnnotationHub::query(ah, c("Homo Sapiens", "release-94"))
#' # snapshotDate(): 2021-09-23
#' # $dataprovider: Ensembl
#' # $species: Homo sapiens
#' # $rdataclass: TwoBitFile, GRanges
#' # additional mcols(): taxonomyid, genome, description, coordinate_1_based,
#' #   maintainer, rdatadateadded, preparerclass, tags,
#' #   rdatapath, sourceurl, sourcetype
#' # retrieve records with, e.g., 'object[["AH64628"]]'
#' #
#' #             title
#' #   AH64628 | Homo_sapiens.GRCh38.94.abinitio.gtf
#' #   AH64629 | Homo_sapiens.GRCh38.94.chr.gtf
#' #   AH64630 | Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf
#' #   AH64631 | Homo_sapiens.GRCh38.94.gtf
#' #   AH65744 | Homo_sapiens.GRCh38.cdna.all.2bit
#' #   AH65745 | Homo_sapiens.GRCh38.dna.primary_assembly.2bit
#' #   AH65746 | Homo_sapiens.GRCh38.dna_rm.primary_assembly.2bit
#' #   AH65747 | Homo_sapiens.GRCh38.dna_sm.primary_assembly.2bit
#' #   AH65748 | Homo_sapiens.GRCh38.ncrna.2bit
#'
#' buildRef(
#'     reference_path = "./Reference_AH",
#'     fasta = "AH65745",
#'     gtf = "AH64631",
#'     genome_type = "hg38"
#' )
#'
#' # Build a SpliceWiz reference, setting chromosome aliases to allow
#' # this reference to process BAM files aligned to UCSC-style genomes:
#'
#' chrom.df <- GenomeInfoDb::genomeStyles()$Homo_sapiens
#'
#' buildRef(
#'     reference_path = "./Reference_UCSC",
#'     fasta = "AH65745",
#'     gtf = "AH64631",
#'     genome_type = "hg38",
#'     chromosome_aliases = chrom.df[, c("Ensembl", "UCSC")]
#' )
#'
#' # One-step generation of SpliceWiz and STAR references, using 4 threads.
#' # NB1: requires a linux-based system with STAR installed.
#' # NB2: A STAR reference genome will be generated in the `STAR` subfolder
#' #      inside the given `reference_path`.
#' # NB3: A custom Mappability Exclusion file will be calculated using STAR
#' #      and will be used to generate the SpliceWiz reference.
#'
#' buildFullRef(
#'     reference_path = "./Reference_with_STAR",
#'     fasta = "genome.fa", gtf = "transcripts.gtf",
#'     genome_type = "",
#'     use_STAR_mappability = TRUE,
#'     n_threads = 4
#' )
#'
#' # NB: the above is equivalent to running the following in sequence:
#'
#' getResources(
#'     reference_path = "./Reference_with_STAR",
#'     fasta = "genome.fa", gtf = "transcripts.gtf"
#' )
#' STAR_buildRef(
#'     reference_path = reference_path,
#'     also_generate_mappability = TRUE,
#'     n_threads = 4
#' )
#' buildRef(
#'     reference_path = "./Reference_with_STAR",
#'     genome_type = ""
#' )
#' }
#' @seealso
#' [Mappability-methods] for methods to calculate low mappability regions\cr\cr
#' [STAR-methods] for a list of STAR wrapper functions\cr\cr
#' \link[AnnotationHub]{AnnotationHub}\cr\cr
#' @name Build-Reference-methods
#' @md
NULL

#' @describeIn Build-Reference-methods Processes / downloads a copy of the genome
#' and gene annotations and stores this in the "resource" subdirectory
#' of the given reference path
#' @export
getResources <- function(
        reference_path = "./Reference",
        fasta = "", gtf = "",
        overwrite = FALSE, force_download = FALSE
) {
    reference_data <- .get_reference_data(
        reference_path = reference_path,
        fasta = fasta, gtf = gtf,
        verbose = TRUE,
        overwrite = overwrite, force_download = force_download,
        pseudo_fetch_fasta = TRUE, pseudo_fetch_gtf = TRUE
    )
}

#' @describeIn Build-Reference-methods First calls \code{getResources()}
#' (if required). Afterwards creates the SpliceWiz reference in the
#' given reference path
#' @export
buildRef <- function(
        reference_path = "./Reference",
        fasta = "", gtf = "", 
        overwrite = FALSE, force_download = FALSE,
        chromosome_aliases = NULL, genome_type = "",
        nonPolyARef = "", MappabilityRef = "", BlacklistRef = "",
        useExtendedTranscripts = TRUE, lowMemoryMode = TRUE
        
    ) {
    .validate_path(reference_path)
    if (!overwrite && 
            file.exists(file.path(reference_path, "SpliceWiz.ref.gz"))) {
        .log("SpliceWiz reference already exists in given directory", "message")
        return()
    }
    extra_files <- .fetch_genome_defaults(reference_path,
        genome_type, nonPolyARef, MappabilityRef, BlacklistRef,
        force_download = force_download)
	
	session <- shiny::getDefaultReactiveDomain()

    N_steps <- 8
    dash_progress("Reading Reference Files", N_steps)
    reference_data <- .get_reference_data(
        reference_path = reference_path,
        fasta = fasta, gtf = gtf, verbose = TRUE,
        overwrite = overwrite, force_download = force_download,
        pseudo_fetch_fasta = lowMemoryMode, pseudo_fetch_gtf = FALSE)

    dash_progress("Processing gtf file", N_steps)
    reference_data$gtf_gr <- .validate_gtf_chromosomes(
        reference_data$genome, reference_data$gtf_gr)
    reference_data$gtf_gr <- .fix_gtf(reference_data$gtf_gr)
    .process_gtf(reference_data$gtf_gr, reference_path)
    extra_files$genome_style <- .gtf_get_genome_style(reference_data$gtf_gr)
    reference_data$gtf_gr <- NULL # To save memory, remove original gtf
    gc()

	# Check here whether fetching from TwoBitFile is problematic
	reference_data$genome <- .check_2bit_performance(reference_path,
		reference_data$genome)		

    dash_progress("Processing introns", N_steps)
    chromosomes <- .convert_chromosomes(chromosome_aliases)
    .process_introns(reference_path, reference_data$genome,
        useExtendedTranscripts)

    dash_progress("Generating SpliceWiz Reference", N_steps)
    .gen_irf(reference_path, extra_files, reference_data$genome, chromosomes)
    gc()

    dash_progress("Annotating IR-NMD", N_steps)
    if(!is.null(session)) {
        shiny::withProgress(message = "Determining NMD Transcripts", {
            .gen_nmd(reference_path, reference_data$genome)
        })
    } else {
        .gen_nmd(reference_path, reference_data$genome)
    }

    dash_progress("Annotating Splice Events", N_steps)
    .gen_splice(reference_path)
    if (file.exists(file.path(reference_path, "fst", "Splice.fst"))) {
        dash_progress("Translating AS Peptides", N_steps)
        .gen_splice_proteins(reference_path, reference_data$genome)
        .log("Splice Annotations finished\n", "message")
    } else {
        dash_progress("No alternate splicing events detected", N_steps)
    }

    message("Reference build finished")
    dash_progress("Reference build finished", N_steps)

    # Prepare a reference-specific cov_data for reference-only plots:
    cov_data <- .prepare_covplot_data(reference_path)
    saveRDS(cov_data, file.path(reference_path, "cov_data.Rds"))

    # Update settings.Rds only after everything is finalised
    settings.list <- readRDS(file.path(reference_path, "settings.Rds"))
    settings.list$genome_type <- genome_type
    settings.list$nonPolyARef <- nonPolyARef
    settings.list$MappabilityRef <- MappabilityRef
    settings.list$BlacklistRef <- BlacklistRef
    settings.list$useExtendedTranscripts <- useExtendedTranscripts
    settings.list$BuildVersion <- buildref_version

    saveRDS(settings.list, file.path(reference_path, "settings.Rds"))
}

#' @describeIn Build-Reference-methods One-step function that fetches resources,
#'   creates a STAR reference (including mappability calculations), then
#'   creates the SpliceWiz reference
#' @export
buildFullRef <- function(
        reference_path,
        fasta, gtf,
        chromosome_aliases = NULL,
        overwrite = FALSE, force_download = FALSE,
        genome_type = genome_type,
        use_STAR_mappability = FALSE,
        nonPolyARef = getNonPolyARef(genome_type),
        BlacklistRef = "",
        useExtendedTranscripts = TRUE,
        n_threads = 4
) {
    if (!overwrite && 
            file.exists(file.path(reference_path, "SpliceWiz.ref.gz"))) {
        .log("SpliceWiz reference already exists in given directory", "message")
        return()
    }

    .validate_STAR_version()

    getResources(reference_path = reference_path,
        fasta = fasta, gtf = gtf,
        overwrite = overwrite, force_download = force_download)

    STAR_buildRef(reference_path = reference_path,
        also_generate_mappability = use_STAR_mappability,
        n_threads = n_threads)

    buildRef(reference_path = reference_path,
        genome_type = genome_type,
        nonPolyARef = nonPolyARef,
        BlacklistRef = BlacklistRef,
        chromosome_aliases = chromosome_aliases,
        useExtendedTranscripts = useExtendedTranscripts)
}


#' @describeIn Build-Reference-methods Returns the path to the BED file containing
#'   coordinates of known non-polyadenylated transcripts for genomes
#'   \code{hg38}, \code{hg19}, \code{mm10} and \code{mm9},
#' @export
getNonPolyARef <- function(genome_type) {
    if (genome_type == "hg38") {
        nonPolyAFile <- system.file(
            "extra-input-files/Human_hg38_nonPolyA_ROI.bed",
            package = "SpliceWiz"
        )
    } else if (genome_type == "hg19") {
        nonPolyAFile <- system.file(
            "extra-input-files/Human_hg19_nonPolyA_ROI.bed",
            package = "SpliceWiz"
        )
    } else if (genome_type == "mm10") {
        nonPolyAFile <- system.file(
            "extra-input-files/Mouse_mm10_nonPolyA_ROI.bed",
            package = "SpliceWiz"
        )
    } else if (genome_type == "mm9") {
        nonPolyAFile <- system.file(
            "extra-input-files/Mouse_mm9_nonPolyA_ROI.bed",
            package = "SpliceWiz"
        )
    } else {
        nonPolyAFile <- ""
    }
    return(nonPolyAFile)
}

################ Functions that may be exported in later releases ##############

Get_Genome <- function(reference_path, validate = TRUE,
        as_DNAStringSet = FALSE) {
    if (validate) .validate_reference(reference_path)
    twobit <- file.path(reference_path, "resource", "genome.2bit")
    if (file.exists(twobit)) {
        genome <- rtracklayer::TwoBitFile(twobit)
    } else if (file.exists(file.path(reference_path, "settings.Rds"))) {
        settings <- readRDS(file.path(reference_path, "settings.Rds"))
        genome <- .fetch_AH(settings$ah_genome, rdataclass = "TwoBitFile")
    } else {
        .log("In Get_Genome, invalid reference_path supplied")
    }
    if (as_DNAStringSet) genome <- rtracklayer::import(genome)
    return(genome)
}

# Returns gzipped GTF file from reference resource folder. For FeatureCounts.
# For STAR, use .STAR_get_GTF instead
Get_GTF_file <- function(reference_path) {
    .validate_reference(reference_path)
    if (file.exists(file.path(reference_path,
            "resource", "transcripts.gtf.gz"))) {
        return(file.path(reference_path, "resource", "transcripts.gtf.gz"))
    } else {
        .log("In Get_GTF_file, invalid reference_path supplied")
    }
}

################################################################################
# Validation functions

.validate_genome_type <- function(genome_type) {
    if (genome_type != "") {
return(TRUE)
}
    .log(paste("In buildRef(),",
        "genome_type not specified.",
        "This should be either one of 'hg38', 'hg19', 'mm10', 'mm9', or",
        "'other'. If 'other', please provide a nonPolyARef file or leave",
        "blank to omit polyA profiling."
    ))
}

.validate_path <- function(reference_path, subdirs = NULL) {
    if ({
        reference_path != "" &&
        tryCatch(
            ifelse(normalizePath(dirname(reference_path)) != "", TRUE, TRUE),
            error = function(e) FALSE
        )
    }) {
        # continue
    } else {
        .log(paste("Error in 'reference_path',",
            paste0("base path of '", reference_path, "' does not exist")
        ))
    }

    base <- normalizePath(dirname(reference_path))
    if (!dir.exists(file.path(base, basename(reference_path))))
        dir.create(file.path(base, basename(reference_path)))

    if (!is.null(subdirs)) {
        for (subdir in subdirs) {
            dir_to_make <- file.path(base, basename(reference_path), subdirs)
            if (!dir.exists(dir_to_make)) dir.create(dir_to_make)
        }
    }
    return(file.path(base, basename(reference_path)))
}

.validate_reference_resource <- function(reference_path, from = "") {
    ref <- normalizePath(reference_path)
    from_str <- ifelse(from == "", "",
        paste("In function", from, ":"))
    if (!dir.exists(ref)) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": this path does not exist"))
    }
    if (!file.exists(file.path(ref, "settings.Rds"))) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": settings.Rds not found"))
    }
    settings.list <- readRDS(file.path(ref, "settings.Rds"))
    if (!("BuildVersion" %in% names(settings.list)) ||
            settings.list[["BuildVersion"]] < buildref_version) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            "SpliceWiz reference is earlier than current version",
            buildref_version))
    }
}

.validate_reference <- function(reference_path, from = "") {
    ref <- normalizePath(reference_path)
    from_str <- ifelse(from == "", "",
        paste("In function", from, ":"))
    if (!dir.exists(ref)) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": this path does not exist"))
    }
    if (!file.exists(file.path(ref, "settings.Rds"))) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": settings.Rds not found"))
    }
    if (!file.exists(file.path(ref, "SpliceWiz.ref.gz"))) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": SpliceWiz.ref.gz not found"))
    }
    settings.list <- readRDS(file.path(ref, "settings.Rds"))
    if (!("BuildVersion" %in% names(settings.list)) ||
            settings.list[["BuildVersion"]] < buildref_version) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            "SpliceWiz reference is earlier than current version",
            buildref_version))
    }
}

.fetch_genome_defaults <- function(reference_path, genome_type,
        nonPolyARef = "", MappabilityRef = "", BlacklistRef = "",
        force_download = FALSE
) {
    if (!is_valid(nonPolyARef)) {
        nonPolyAFile <- getNonPolyARef(genome_type)
        nonPolyAFile <- .parse_valid_file(nonPolyAFile,
            "Reference generated without non-polyA reference")
    } else {
        nonPolyAFile <- .parse_valid_file(nonPolyARef,
            "Reference generated without non-polyA reference",
            force_download = force_download)
    }
    map_path <- file.path(normalizePath(reference_path), "Mappability")
    map_file <- file.path(map_path, "MappabilityExclusion.bed.gz")
    if (is_valid(MappabilityRef)) {
        MappabilityFile <- .parse_valid_file(MappabilityRef,
            force_download = force_download)
    } else if (file.exists(map_file)) {
        MappabilityFile <- .parse_valid_file(map_file)
    } else if (genome_type %in% c("hg38", "hg19", "mm9", "mm10")) {
        map.gz <- get_mappability_exclusion(
            genome_type, as_type = "bed.gz", path = map_path, overwrite = TRUE)
        if(!is.null(map.gz)) {
            MappabilityFile <- .parse_valid_file(map.gz)
        } else {
            MappabilityFile <- ""
        }
    } else {
        MappabilityFile <- .parse_valid_file(MappabilityRef,
            "Reference generated without Mappability reference")
    }
    BlacklistFile <-
        .parse_valid_file(BlacklistRef,
            "Reference generated without Blacklist exclusion",
            force_download = force_download)

    # Check files are valid BED files; fail early if not
    .check_is_BED(nonPolyAFile)
    .check_is_BED(MappabilityFile)
    .check_is_BED(BlacklistFile)

    final <- list(
        nonPolyAFile = nonPolyAFile, MappabilityFile = MappabilityFile,
        BlacklistFile = BlacklistFile, genome_type = genome_type
    )
    return(final)
}

.check_is_BED <- function(filename) {
    if (is_valid(filename)) {
        tryCatch(rtracklayer::import.bed(filename, "bed"),
            error = function(x) {
                .log(paste(filename, "is not a BED file"))
        })
    }
    return()
}

.convert_chromosomes <- function(chromosome_aliases) {
    if (is.null(chromosome_aliases)) return(NULL)
    df <- as.data.frame(chromosome_aliases)
    df <- df[!duplicated(df[, 1]), ]
    df <- df[!duplicated(df[, 2]), ]
    df <- df[as.vector(df[, 1]) != "", ]
    df <- df[as.vector(df[, 2]) != "", ]
    df <- as.data.table(df)
    colnames(df) <- c("Original", "New")
    return(df)
}

################################################################################
# Sub

.get_reference_data <- function(reference_path, fasta, gtf,
        verbose = TRUE, overwrite = FALSE, force_download = FALSE,
        pseudo_fetch_fasta = FALSE, pseudo_fetch_gtf = FALSE
) {

    # Checks fasta or gtf files exist if omitted, or are valid URLs
    .validate_path(reference_path, subdirs = "resource")
    if (!is_valid(fasta)) {
        twobit <- file.path(reference_path, "resource", "genome.2bit")
        if (!file.exists(twobit)) .log(paste(fasta, "doesn't exist"))
    }
    if (!is_valid(gtf)) {
        gtf <- file.path(reference_path, "resource", "transcripts.gtf.gz")
        if (!file.exists(gtf)) .log(paste(gtf, "doesn't exist"))
    }
    # URLS are checked at BiocFileCache step contained within .parse_valid_file

    fasta_use <- gtf_use <- ah_genome_use <- ah_gtf_use <- ""
    if (.is_AH_pattern(fasta)) {
        ah_genome_use <- fasta
    } else {
        fasta_use <- fasta
    }
    if (.is_AH_pattern(gtf)) {
        ah_gtf_use <- gtf
    } else {
        gtf_use <- gtf
    }
	session <- shiny::getDefaultReactiveDomain()
	fetch_fasta_args <- list(
        reference_path = reference_path,
        fasta = fasta_use, ah_genome = ah_genome_use,
        verbose = verbose, overwrite = overwrite,
        force_download = force_download, pseudo_fetch = pseudo_fetch_fasta
    )
	fetch_gtf_args <- list(
        gtf = gtf_use, ah_transcriptome = ah_gtf_use,
        reference_path = reference_path,
        verbose = verbose, overwrite = overwrite,
        force_download = force_download, pseudo_fetch = pseudo_fetch_gtf
    )
	genome <- gtf_gr <- NULL
	
	# wrap nested progress bars for shiny
    if(!is.null(session)) {
        shiny::withProgress(message = "Loading genome", {
            genome <- do.call(.fetch_fasta, fetch_fasta_args)
        })
        shiny::withProgress(message = "Loading gene annotations", {
            gtf_gr <- do.call(.fetch_gtf, fetch_gtf_args)
        })
    } else {
		genome <- do.call(.fetch_fasta, fetch_fasta_args)
		gtf_gr <- do.call(.fetch_gtf, fetch_gtf_args)
    }

    # Save Resource details to settings.Rds:
    settings.list <- list(fasta_file = fasta_use, gtf_file = gtf_use,
        ah_genome = ah_genome_use, ah_transcriptome = ah_gtf_use,
        reference_path = reference_path
    )
    settings.list$BuildVersion <- buildref_version
    saveRDS(settings.list, file.path(reference_path, "settings.Rds"))

    settings.list <- readRDS(file.path(reference_path, "settings.Rds"))
    final <- list(genome = genome, gtf_gr = gtf_gr)
    return(final)
}

.is_AH_pattern <- function(word) {
    if (substr(word, 1, 2) == "AH" && !file.exists(word)) return(TRUE)
    return(FALSE)
}


################################################################################

.fetch_genome_as_required <- function(genome, pseudo_fetch) {
    if (!pseudo_fetch) {
        .log("Importing genome sequences to memory...", "message",
            appendLF = FALSE)
        genome <- import(genome) # Fetch as DNAStringSet - avoid TwoBit lag
        message("done")
    }
    return(genome)
}

# If fasta is a web resource
#   - downloads this and rename as genome.fa inside "resource" dir
#   - make a genome.2bit file
# If fasta is a file
#   - make a copy and rename as genome.fa inside "resource" dir
#   - make a genome.2bit file
# If ah_genome is not empty:
#   - fetch AnnotationHub reference
#   - create a local genome.2bit file for portability
.fetch_fasta <- function(
        reference_path = "./Reference",
        fasta = "", ah_genome = "",
        verbose = TRUE, overwrite = FALSE, force_download = FALSE,
        pseudo_fetch = FALSE
) {
    if (ah_genome != "") {
		N_steps <- 3
		dash_progress("Retrieving AnnotationHub resource", N_steps)
        genome <- .fetch_fasta_ah(ah_genome, verbose = verbose)
		dash_progress("Saving copy of genome as TwoBit file", N_steps)
        .fetch_fasta_save_2bit(genome, reference_path, overwrite)
		dash_progress("Loading genome into memory", N_steps)
        genome <- .fetch_genome_as_required(genome, pseudo_fetch)
        return(genome)
    } else if (fasta == "") {
        twobit <- file.path(reference_path, "resource", "genome.2bit")
        if (file.exists(twobit)) {
			N_steps <- 1
			dash_progress("Loading genome into memory", N_steps)
            .log("Connecting to genome TwoBitFile...", "message",
                appendLF = FALSE)
            genome <- Get_Genome(reference_path, validate = FALSE,
                as_DNAStringSet = !pseudo_fetch)
            message("done")
            return(genome)
        } else {
            .log("Resource genome is not available; `fasta` parameter required")
        }
    } else {
        # If no overwrite, quickly return genome.2bit if exists
        if (!overwrite) {
            twobit <- file.path(reference_path, "resource", "genome.2bit")
            if (file.exists(twobit)) {
				N_steps <- 1
				dash_progress("Loading genome into memory", N_steps)
                .log("Connecting to genome TwoBitFile...", "message",
                    appendLF = FALSE)
                genome <- Get_Genome(reference_path, validate = FALSE,
                    as_DNAStringSet = !pseudo_fetch)
                message("done")
                return(genome)
            }
        }
		N_steps <- 4
        .log("Converting FASTA to local TwoBitFile...", "message",
            appendLF = FALSE)
		dash_progress("Downloading genome, if required...", N_steps)
        fasta_file <- .parse_valid_file(fasta, force_download = force_download)
        if (!file.exists(fasta_file))
            .log(paste("In .fetch_fasta(),",
                "Given genome fasta file", fasta, "not found"))

		dash_progress("Loading genome into memory", N_steps)
        genome <- .fetch_fasta_file(fasta_file)
		dash_progress("Saving copy of genome as TwoBit file", N_steps)
        .fetch_fasta_save_2bit(genome, reference_path, overwrite)
        message("done")

		rm(genome)
		gc() # free memory before re-import
		.log("Connecting to genome TwoBitFile...", "message", appendLF = FALSE)
		dash_progress("Reloading genome from TwoBit file...", N_steps)
		genome <- Get_Genome(reference_path, validate = FALSE,
			as_DNAStringSet = !pseudo_fetch)
		# TwoBitFile's getSeq is slow on some linux systems (don't know why)
		# Importing TwoBitFile as a proper DNAStringSet
		message("done")
        return(genome)
    }
}

# Fetch the AnnotationHub resource and return as a genome object
.fetch_fasta_ah <- function(ah_genome, verbose = TRUE) {
    if (substr(ah_genome, 1, 2) != "AH")
        .log("Given genome AnnotationHub reference is incorrect")

    genome <- .fetch_AH(ah_genome, verbose = verbose, rdataclass = "TwoBitFile",
        as_DNAStringSet = FALSE)
}

.fetch_fasta_file <- function(fasta_file) {
    # .log("Importing genome into memory...", "message", appendLF = FALSE)
    genome <- Biostrings::readDNAStringSet(fasta_file)
    # message("done")
    return(genome)
}

.fetch_fasta_save_fasta <- function(genome, reference_path, overwrite) {
    genome.fa <- file.path(reference_path, "resource", "genome.fa")
    if (overwrite || !file.exists(genome.fa)) {
        .log("Saving local copy as FASTA...", "message", appendLF = FALSE)
        if (overwrite && file.exists(genome.fa)) {
            file.remove(genome.fa)
        }
        rtracklayer::export(genome, genome.fa, "fasta")
        message("done")
    }
}

.fetch_fasta_save_2bit <- function(genome, reference_path, overwrite) {
    genome.2bit <- file.path(reference_path, "resource", "genome.2bit")
    if (is(genome, "TwoBitFile") && file.exists(genome.2bit) &&
            normalizePath(rtracklayer::path(genome)) ==
            normalizePath(genome.2bit)) {
        return()
    } # prevent self-writing
    if (overwrite || !file.exists(genome.2bit)) {
        # .log("Saving genome as TwoBitFile...", "message", appendLF = FALSE)
        if (overwrite && file.exists(genome.2bit)) {
            file.remove(genome.2bit)
        }
        if (is(genome, "TwoBitFile") && 
                file.exists(rtracklayer::path(genome))) {
            file.copy(rtracklayer::path(genome), genome.2bit)
        } else {
            rtracklayer::export(genome, genome.2bit, "2bit")
        }
        # message("done")
    }
}

################################################################################

.fetch_gtf <- function(
        reference_path = "./Reference",
        gtf = "", ah_transcriptome = "",
        verbose = TRUE, overwrite = FALSE, force_download = FALSE,
        pseudo_fetch = FALSE
) {
    r_path <- file.path(reference_path, "resource")
    gtf_path <- file.path(r_path, "transcripts.gtf.gz")
    if (ah_transcriptome != "") {
		N_steps <- 1
		dash_progress("Retrieving AnnotationHub resource", N_steps)
        gtf_gr <- .fetch_AH(ah_transcriptome, verbose = verbose,
            pseudo_fetch = pseudo_fetch)
        if (overwrite || !file.exists(gtf_path)) {
            cache_loc <- .fetch_AH_cache_loc(ah_transcriptome,
                rdataclass = "GRanges", verbose = verbose)
            if (file.exists(cache_loc)) {
                if (file.exists(gtf_path)) file.remove(gtf_path)
                file.copy(cache_loc, gtf_path)
            } # Copy file from cache if exists
        }
        return(gtf_gr)
    } else if (gtf == "") {
        if (file.exists(gtf_path)) {
			N_steps <- 1
			dash_progress("Loading gene annotation into memory", N_steps)
            .log("Reading source GTF file...", "message", appendLF = FALSE)
            gtf_gr <- rtracklayer::import(gtf_path, "gtf")
            message("done")
			return(gtf_gr)
        } else {
            .log(paste("Resource gene annotation is not available;",
				"`gtf` parameter required"))
        }
		
    } else {
        # If no overwrite, quickly return GTF if exists
        if (!overwrite) {
            if (file.exists(gtf_path)) {
				N_steps <- 1
				dash_progress("Loading gene annotation into memory", N_steps)
				.log("Reading source GTF file...", "message", appendLF = FALSE)
				gtf_gr <- rtracklayer::import(gtf_path, "gtf")
				message("done")
				return(gtf_gr)
            }
        }
		
		N_steps <- 3
		dash_progress("Downloading gene annotations, if required...", N_steps)
        gtf_file <- .parse_valid_file(gtf, force_download = force_download)
        if (!file.exists(gtf_file)) {
            .log(paste("In .fetch_gtf(),",
                "Given transcriptome gtf file", gtf, "not found"))
        }
		dash_progress("Making local copy, if required...", N_steps)
        if (!file.exists(gtf_path) ||
                normalizePath(gtf_file) != normalizePath(gtf_path)) {
            if (overwrite || !file.exists(gtf_path)) {
                .log("Making local copy of GTF file...", "message",
                    appendLF = FALSE)
                if (substr(gtf_file, nchar(gtf_file) - 2,
                        nchar(gtf_file)) == ".gz") {
                    if (file.exists(gtf_path)) file.remove(gtf_path)
                    file.copy(gtf_file, gtf_path)
                } else {
                    gzip(filename = gtf_file, destname = gtf_path,
                        remove = FALSE)
                }
                message("done")
            }
        }
		dash_progress("Loading annotations into memory", N_steps)
        if (!pseudo_fetch) {
            .log("Reading source GTF file...", "message", appendLF = FALSE)
            gtf_gr <- rtracklayer::import(gtf_file, "gtf")
            message("done")
        } else {
            gtf_gr <- NULL
        }
        return(gtf_gr)
    }
}

# Check some chromosomes exist between gtf and genome
.validate_gtf_chromosomes <- function(genome, gtf_gr) {
    chrOrder <- names(seqinfo(genome))
    if (!any(as.character(GenomicRanges::seqnames(gtf_gr)) %in% chrOrder)) {
        .log(paste("In .validate_gtf_chromosomes(),",
            "Chromosomes in genome and gene annotation does not match",
            "likely incompatible FASTA and GTF file"))
    }
    seqlevels(gtf_gr, pruning.mode = "tidy") <- chrOrder
    return(gtf_gr)
}

################################################################################

.fetch_AH_cache_loc <- function(ah_record_name,
    rdataclass = c("GRanges", "TwoBitFile"),
    localHub = FALSE, ah = AnnotationHub(localHub = localHub),
    verbose = FALSE
) {
    rdataclass <- match.arg(rdataclass)
    if (!substr(ah_record_name, 1, 2) == "AH")
        .log(paste(ah_record_name,
            "does not appear to be a valid AnnotationHub record name"))

    if (!(ah_record_name %in% names(ah)))
        .log(paste(ah_record_name,
            "is not found in AnnotationHub index.",
            "Perhaps check online connection or record name"))

    ah_record <- ah[names(ah) == ah_record_name]
    if (ah_record$rdataclass != rdataclass)
        .log(paste(ah_record_name,
            "is of type", ah_record$rdataclass,
            "and not of expected:", rdataclass))

    if (verbose)
        .log(paste("Downloading", rdataclass,
            "from AnnotationHub, if required..."),
            "message", appendLF = FALSE)

    cache_loc <- AnnotationHub::cache(ah_record)
    if (verbose) message("done")
    if (!file.exists(cache_loc))
        .log("AnnotationHub cache error - asset not found")

    return(cache_loc)
}

.fetch_AH <- function(ah_record_name, rdataclass = c("GRanges", "TwoBitFile"),
        localHub = FALSE, ah = AnnotationHub(localHub = localHub),
        as_DNAStringSet = FALSE, verbose = FALSE,
        pseudo_fetch = pseudo_fetch
) {
    rdataclass <- match.arg(rdataclass)
    cache_loc <- .fetch_AH_cache_loc(ah_record_name, rdataclass,
        localHub, ah, verbose)
    if (rdataclass == "GRanges") {
        if (!pseudo_fetch) {
            if (verbose) .log("Importing to memory as GRanges object...",
                "message", appendLF = FALSE)

            gtf <- rtracklayer::import(cache_loc, "gtf")
            if (verbose) message("done")
            return(gtf)
        }
        return(NULL)
    } else if (rdataclass == "TwoBitFile") {
        if (verbose) .log("Importing to memory as TwoBitFile object...",
            "message", appendLF = FALSE)

        twobit <- rtracklayer::TwoBitFile(cache_loc)
        if (verbose) message("done")
        if (as_DNAStringSet) {
            .log("Importing genome into memory...", "message", appendLF = FALSE)
            genome <- rtracklayer::import(twobit)
            message("done")
            return(genome)
        }
        return(twobit)
    }
}

# debug function
.get_SpliceWiz_cache <- function() {
    cache <- tools::R_user_dir(package = "SpliceWiz", which = "cache")
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    bfc
}

.parse_valid_file <- function(file, msg = "", force_download = FALSE) {
    if (!is_valid(file)) {
        .log(msg, type = "message")
        return("")
    } else if (any(startsWith(file, c("http", "ftp")))) {
        url <- file
        # TODO: test URLs here
        cache <- tools::R_user_dir(package = "SpliceWiz", which = "cache")
        bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
        res <- BiocFileCache::bfcquery(bfc, url, "fpath", exact = TRUE)
        if (nrow(res) > 0 & !force_download)
            return(paste(cache, res$rpath[nrow(res)], sep = "/"))

        # either force_download == TRUE or nrow(res) == 0
        path <- tryCatch(BiocFileCache::bfcadd(bfc, url),
            error = function(err) {
                .log(paste("Web resource not accessible -", url), "message")
                return(NA)
            }
        )
        if (identical(path, NA)) {
            # fetch local copy if available
            if (nrow(res) == 0) return("")
            .log("Returning local copy from cache", "message")
            return(paste(cache, res$rpath[nrow(res)], sep = "/"))
        }
        # remove prior versions from cache to remove glut
        res <- BiocFileCache::bfcquery(bfc, url, "fpath", exact = TRUE)
        if(nrow(res) > 1)
            BiocFileCache::bfcremove(bfc, res$rid[-nrow(res)])
        return(paste(cache, res$rpath[nrow(res)], sep = "/"))
    } else if (!file.exists(file)) {
        .log(paste(file, "not found.", msg), type = "message")
        return("")
    } else if (file.exists(file)) {
        return(file)
    } else {
        .log(msg, type = "message")
        return("")
    }
}

# Various GTF fixes
# - fix gene / transcript names with '/' (which breaks processBAM C++ code)
# - also fix missing gene_name and transcript_names in newer Ensembl refs
.fix_gtf <- function(gtf_gr) {

    if ("gene_name" %in% names(S4Vectors::mcols(gtf_gr))) {
        gtf_gr$gene_name[is.na(gtf_gr$gene_name)] <-
            gtf_gr$gene_id[is.na(gtf_gr$gene_name)]
        gtf_gr$gene_name <- gsub("/", "_", gtf_gr$gene_name)
    } else {
        gtf_gr$gene_name <- gtf_gr$gene_id
    }

    if ("transcript_name" %in% names(S4Vectors::mcols(gtf_gr))) {
        gtf_gr$transcript_name[is.na(gtf_gr$transcript_name)] <-
            gtf_gr$transcript_id[is.na(gtf_gr$transcript_name)]
        gtf_gr$transcript_name <- gsub("/", "_", gtf_gr$transcript_name)
    } else {
        gtf_gr$transcript_name <- gtf_gr$transcript_id
    }
    if (!("gene_biotype" %in% names(S4Vectors::mcols(gtf_gr)))) {
        gtf_gr$gene_biotype <- "protein_coding"
    }
    if (!("transcript_biotype" %in% names(S4Vectors::mcols(gtf_gr)))) {
        gtf_gr$transcript_biotype <- "protein_coding"
    }
    if (!("transcript_support_level" %in% names(S4Vectors::mcols(gtf_gr)))) {
        gtf_gr$transcript_support_level <- 1
    }

    return(gtf_gr)
}

.gtf_get_genome_style <- function(gtf_gr) {
    seqnames <- names(seqinfo(gtf_gr))
    UCSC <- any(seqnames %in% genomeStyles("Homo sapiens")$UCSC)
    Ensembl <- any(seqnames %in% genomeStyles("Homo sapiens")$Ensembl)
    if (UCSC == Ensembl) {
        return("")
    } else if (UCSC) {
        return("UCSC")
    } else {
        return("Ensembl")
    }
}

################################################################################
# Sub

.process_gtf <- function(gtf_gr, reference_path) {
    # Create "fst" subdirectory if not exists
    .validate_path(reference_path, subdirs = "fst")
    .log("Processing gtf file...", "message")
    message("...genes")
    Genes_group <- .process_gtf_genes(gtf_gr, reference_path)
    message("...transcripts")
    .process_gtf_transcripts(gtf_gr, reference_path)
    message("...CDS")
    .process_gtf_misc(gtf_gr, reference_path)
    message("...exons")
    .process_gtf_exons(gtf_gr, reference_path, Genes_group)
    message("done")
}

# Processes Genes
# - includes computation of gene groups
.process_gtf_genes <- function(gtf_gr, reference_path) {
    Genes <- gtf_gr[gtf_gr$type == "gene"]

    # If genes are not annotated by "type" column, then have to do it manually
    if (length(Genes) == 0) {
        gene_cols <- c("seqnames", "strand",
            "gene_id", "gene_name", "gene_biotype")
        Genes <- as.data.table(gtf_gr)
        Genes <- Genes[, c("start", "end", "width") := list(
            min(get("start")), max(get("end")),
            max(get("end")) - min(get("start")) + 1
        ), by = gene_cols]
        Genes <- unique(Genes, by = gene_cols)
        Genes$type <- "gene"
        Genes <- .grDT(Genes, keep.extra.columns = TRUE)
        if (length(Genes) == 0) .log("No genes detected in reference!")
    }

    Genes <- GenomeInfoDb::sortSeqlevels(Genes)
    Genes <- sort(Genes)
    Genes$gene_display_name <- paste0(Genes$gene_name, " (", Genes$gene_id, ")")

    # Annotate gene_groups_stranded
    Genes_group.stranded <- as.data.table(reduce(Genes))
    setorder(Genes_group.stranded, seqnames, start, strand)

    Genes_group.stranded[, c("gene_group") := .I]
    OL <- findOverlaps(
        Genes, .grDT(Genes_group.stranded)
    )
    Genes$gene_group_stranded[from(OL)] <-
        Genes_group.stranded$gene_group[to(OL)]

    # Annotate gene_groups_unstranded
    Genes_group.unstranded <- as.data.table(reduce(Genes, ignore.strand = TRUE))
    setorder(Genes_group.unstranded, seqnames, start)
    Genes_group.unstranded[, c("gene_group") := .I]
    OL <- findOverlaps(
        Genes,
        .grDT(Genes_group.unstranded, ignore.strand = TRUE)
    )
    Genes$gene_group_unstranded[from(OL)] <-
        Genes_group.unstranded$gene_group[to(OL)]

    write.fst(as.data.frame(Genes),
        file.path(reference_path, "fst", "Genes.fst")
    )
    final <- list(
        stranded = Genes_group.stranded,
        unstranded = Genes_group.unstranded
    )
    return(final)
}

.process_gtf_transcripts <- function(gtf_gr, reference_path) {
    Transcripts <- gtf_gr[gtf_gr$type == "transcript"]

    # If transcript are not annotated by "type" column, then do manually
    if (length(Transcripts) == 0) {
        tx_cols <- c("seqnames", "strand",
            "gene_id", "gene_name", "gene_biotype",
            "transcript_id", "transcript_name", "transcript_biotype")
        Transcripts <- as.data.table(gtf_gr)
        Transcripts <- Transcripts[, c("start", "end", "width") := list(
            min(get("start")), max(get("end")),
            max(get("end")) - min(get("start")) + 1
        ), by = tx_cols]
        Transcripts <- unique(Transcripts, by = tx_cols)
        Transcripts$type <- "transcript"
        Transcripts <- .grDT(Transcripts, keep.extra.columns = TRUE)
        if (length(Transcripts) == 0) 
            .log("No transcripts detected in reference!")
    }

    Transcripts <- GenomeInfoDb::sortSeqlevels(Transcripts)
    Transcripts <- sort(Transcripts)

    # Fix gene_biotype and transcript_biotype tags
    if ("gene_biotype" %in% names(mcols(Transcripts))) {
        # do nothing
    } else if ("gene_type" %in% names(mcols(Transcripts))) {
        colnames(mcols(Transcripts))[which(colnames(mcols(Transcripts)) ==
            "gene_type")] <- "gene_biotype"
    } else {
        mcols(Transcripts)$gene_biotype <- "protein_coding"
    }
    if ("transcript_biotype" %in% names(mcols(Transcripts))) {
        # do nothing
    } else if ("transcript_type" %in% names(mcols(Transcripts))) {
        colnames(mcols(Transcripts))[which(colnames(mcols(Transcripts)) ==
            "transcript_type")] <- "transcript_biotype"
    } else {
        mcols(Transcripts)$transcript_biotype <- "protein_coding"
    }
    if ("transcript_support_level" %in% names(mcols(Transcripts))) {
        Transcripts$transcript_support_level <-
            tstrsplit(Transcripts$transcript_support_level, split = " ")[[1]]
        Transcripts$transcript_support_level[
            is.na(Transcripts$transcript_support_level)
        ] <- "NA"
    }
    write.fst(as.data.frame(Transcripts),
        file.path(reference_path, "fst", "Transcripts.fst")
    )
}

.process_gtf_misc <- function(gtf_gr, reference_path) {
    # Proteins
    Proteins <- gtf_gr[gtf_gr$type == "CDS"]
    if (length(Proteins) == 0) {
        .log("No CDS (proteins) detected in reference!")
    } # Is this critical to SpliceWiz function?

    Proteins <- GenomeInfoDb::sortSeqlevels(Proteins)
    Proteins <- sort(Proteins)
    write.fst(
        as.data.frame(Proteins),
        file.path(reference_path, "fst", "Proteins.fst")
    )
    # Misc
    gtf.misc <- gtf_gr[!gtf_gr$type %in% c("gene", "transcript", "exon", "CDS")]
    if (length(gtf.misc) == 0) {
        .log("No start / stop codons detected in reference!")
    }
    gtf.misc <- GenomeInfoDb::sortSeqlevels(gtf.misc)
    gtf.misc <- sort(gtf.misc)
    write.fst(
        as.data.frame(gtf.misc),
        file.path(reference_path, "fst", "Misc.fst")
    )
}

.process_gtf_exons <- function(gtf_gr, reference_path, Genes_group) {
    Exons <- gtf_gr[gtf_gr$type == "exon"]
    if (length(Exons) == 0) .log("No exons detected in reference!")

    Exons <- GenomeInfoDb::sortSeqlevels(Exons)
    Exons <- sort(Exons)

    # transcript_biotype is very important field.
    #   If Gencode, this is transcript_type.
    #   In rare case we do not have this field
    #   This next bit ensures transcript_biotype exists.
    if ("transcript_biotype" %in% names(mcols(Exons))) {
    } else if ("transcript_type" %in% names(mcols(Exons))) {
        colnames(mcols(Exons))[
            which(colnames(mcols(Exons)) == "transcript_type")
        ] <- "transcript_biotype"
    } else {
        mcols(Exons)$transcript_biotype <- "protein_coding"
    }

    # Assign gene groups then bake exon-groups into Exons
    tmp.Exons_group.stranded <- .process_exon_groups(
        Exons, Genes_group, stranded = TRUE)
    tmp.Exons_group.unstranded <- .process_exon_groups(
        Exons, Genes_group, stranded = FALSE)

    # Now annotate all exons in Exons with the gene and exon groups
    OL <- findOverlaps(Exons, .grDT(tmp.Exons_group.stranded))
    Exons$gene_group_stranded[from(OL)] <-
        tmp.Exons_group.stranded$gene_group[to(OL)]
    Exons$exon_group_stranded[from(OL)] <-
        tmp.Exons_group.stranded$exon_group[to(OL)]

    OL <- findOverlaps(Exons, .grDT(tmp.Exons_group.unstranded,
            ignore.strand = TRUE
        )
    )
    Exons$gene_group_unstranded[from(OL)] <-
        tmp.Exons_group.unstranded$gene_group[to(OL)]
    Exons$exon_group_unstranded[from(OL)] <-
        tmp.Exons_group.unstranded$exon_group[to(OL)]

    write.fst(as.data.frame(Exons),
        file.path(reference_path, "fst", "Exons.fst"))
    write.fst(
        rbind(tmp.Exons_group.stranded, tmp.Exons_group.unstranded),
        file.path(reference_path, "fst", "Exons.Group.fst")
    )
}

.process_exon_groups <- function(Exons, Genes_group, stranded = TRUE) {

    # Annotated IR transcripts exons are not regarded as exons
    tmp.exons.exclude <- Exons[!grepl("intron", Exons$transcript_biotype)]
    tmp.Exons_group <- as.data.table(reduce(tmp.exons.exclude,
        ignore.strand = !stranded
    ))
    GG <- Genes_group[[ifelse(stranded, "stranded", "unstranded")]]

    OL <- findOverlaps(
        .grDT(tmp.Exons_group),
        .grDT(GG, ignore.strand = !stranded)
    )
    tmp.Exons_group$gene_group[from(OL)] <- GG$gene_group[to(OL)]

    # Some retained_intron transcripts have terminal exons lying outside of
    #   main transcripts. Include these also
    tmp.exons.exclude.span <- split(
        .grDT(tmp.Exons_group),
        tmp.Exons_group$gene_group
    )
    tmp.exons.exclude.span <-
        unlist(range(tmp.exons.exclude.span), use.names = TRUE)
    tmp.exons.RI <- Exons[grepl("intron", Exons$transcript_biotype)]
    if (length(tmp.exons.RI) > 0) {
        OL <- findOverlaps(
            tmp.exons.RI,
            tmp.exons.exclude.span,
            ignore.strand = !stranded
        )
        tmp.exons.RI <- tmp.exons.RI[-from(OL)]
        tmp.exons.exclude <- c(tmp.exons.exclude, tmp.exons.RI)

        tmp.Exons_group <- as.data.table(reduce(tmp.exons.exclude,
            ignore.strand = !stranded
        ))
        OL <- findOverlaps(
            .grDT(tmp.Exons_group),
            .grDT(GG, ignore.strand = !stranded)
        )
        tmp.Exons_group$gene_group[from(OL)] <-
            GG$gene_group[to(OL)]
    }
    setorder(tmp.Exons_group, seqnames, start, strand)
    tmp.Exons_group[, c("exon_group") :=
        data.table::rowid(get("gene_group"))]
    if (stranded) {
        tmp.Exons_group[get("strand") == "-",
            c("exon_group") := max(get("exon_group")) + 1 - get("exon_group"),
            by = "gene_group"
        ]
    }
    return(tmp.Exons_group)
}

# Fetch first and last 1000 exon sequences, then evaluates the ratio of times
#   taken to fetch. A ratio of > 3 is evidence there is position-dependent
#   loading lag (which is observed on some systems)
.check_2bit_performance <- function(reference_path, genome) {
	if(is(genome, "TwoBitFile")) {
		Exons <- as.data.table(
			read.fst(file.path(reference_path, "fst", "Exons.fst")),
		)
		if(nrow(Exons) > 1000) {
			gr_start <- .grDT(Exons[seq_len(1000)])
			gr_end <- .grDT(Exons[seq(nrow(Exons) - 999, nrow(Exons))])
			bench_start <- system.time(getSeq(genome, gr_start))
			bench_end <- system.time(getSeq(genome, gr_end))
			# message("TwoBit fetch benchmark start: ", unname(bench_start[3]), 
				# ", end: ", unname(bench_end[3]))
			if(bench_start[3] > 0 && bench_end[3] / bench_start[3] > 3) {
				.log(paste("SpliceWiz detected inefficient TwoBit retrieval,",
					" importing genome as DNAStringSet"),
					"message")
				return(rtracklayer::import(genome))
			}
		}
	}
	return(genome)
}

################################################################################
# Sub

.process_introns <- function(reference_path, genome,
        useExtendedTranscripts = TRUE) {
    .log("Processing introns...", "message")

    message("...data")
    data <- .process_introns_data(reference_path, genome, 
        useExtendedTranscripts)
    gc()
    data[["candidate.introns"]] <- .process_introns_annotate(
        data[["candidate.introns"]], data[["Transcripts"]], genome,
        data[["Proteins"]], data[["Exons"]]
    )
    message("...defining flanking exon islands")
    data[["candidate.introns"]] <- .process_introns_group(
        data[["candidate.introns"]], data[["Exons_group.stranded"]],
        data[["Exons_group.unstranded"]]
    )
    gc()
    write.fst(data[["candidate.introns"]],
        file.path(reference_path, "fst", "junctions.fst"))
    message("done")
}

# Import data for intron processing; create list of candidate.introns
.process_introns_data <- function(reference_path, genome,
        useExtendedTranscripts = TRUE) {
    Exons <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.fst")),
    )
    Transcripts <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Transcripts.fst")),
    )
    Proteins <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Proteins.fst")),
    )
    Exons_group <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.Group.fst")),
    )
    Exons_group.stranded <- Exons_group[get("strand") != "*"]
    Exons_group.unstranded <- Exons_group[get("strand") == "*"]

    if (useExtendedTranscripts == FALSE) {
        candidate.transcripts <- Exons[get("transcript_biotype") %in%
            c("processed_transcript", "protein_coding")]
    } else {
        candidate.transcripts <- Exons
    }
    candidate.introns <- as.data.table(
        .grlGaps(
            split(.grDT(candidate.transcripts),
            candidate.transcripts$transcript_id)
        )
    )
    candidate.introns[, c("group") := NULL]
    setnames(candidate.introns, "group_name", "transcript_id")
    setorderv(candidate.introns, c("seqnames", "start", "end", "strand"))

    final <- list(
        Exons = Exons, Transcripts = Transcripts, Proteins = Proteins,
        Exons_group.stranded = Exons_group.stranded,
        Exons_group.unstranded = Exons_group.unstranded,
        candidate.introns = candidate.introns
    )
    return(final)
}

#############################################################################

# Annotate particulars for given junctions / introns
.process_introns_annotate <- function(candidate.introns,
        Transcripts, genome, Proteins, Exons) {
    # Annotating Intron number, gene/transcript names, biotype:
    message("...basic annotations")
    candidate.introns <- .process_introns_annotate_basics(
        candidate.introns, Transcripts)

    # Grab splice motifs
    message("...splice motifs")
    candidate.introns <- .process_introns_annotate_splice_motifs(
        candidate.introns, genome)

    # Annotate TSL, protein_id, ccds_id:
    message("...other info")
    candidate.introns <- .process_introns_annotate_others(
        candidate.introns, Transcripts, Proteins, Exons
    )
    gc()
    return(candidate.introns)
}

.process_introns_annotate_basics <- function(
        candidate.introns, Transcripts) {
    # Intron number
    candidate.introns[, c("intron_number") :=
        data.table::rowid(get("transcript_id"))]
    candidate.introns[get("strand") == "-",
        c("intron_number") :=
            max(get("intron_number")) + 1 - get("intron_number"),
        by = "transcript_id"
    ]

    # Name introns by number
    candidate.introns[, c("intron_id") :=
        paste0(get("transcript_id"), "_Intron", get("intron_number"))]

    # Annotate introns with gene and transcript names and biotype
    candidate.introns[Transcripts,
        on = "transcript_id",
        c("gene_name", "gene_id", "transcript_name", "transcript_biotype") :=
            list(
                get("i.gene_name"), get("i.gene_id"),
                get("i.transcript_name"), get("i.transcript_biotype")
            )
    ]
    return(candidate.introns)
}

.process_introns_annotate_splice_motifs <- function(
        candidate.introns, genome) {

    donor.introns <- data.frame(
        seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+",
            candidate.introns$start, candidate.introns$end - 1
        ),
        stop = ifelse(candidate.introns$strand == "+",
            candidate.introns$start + 1, candidate.introns$end
        ),
        strand = candidate.introns$strand
    )
    donor_seq <- getSeq(genome, .grDT(donor.introns))
    acceptor.introns <- data.frame(
        seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+",
            candidate.introns$end - 1, candidate.introns$start
        ),
        stop = ifelse(candidate.introns$strand == "+",
            candidate.introns$end, candidate.introns$start + 1
        ),
        strand = candidate.introns$strand
    )
    acceptor_seq <- getSeq(genome, .grDT(acceptor.introns))
    candidate.introns$splice_motif <- paste0(donor_seq, acceptor_seq)

    return(candidate.introns)
}

.process_introns_annotate_others <- function(
        candidate.introns, Transcripts, Proteins, Exons) {

    # Annotate TSL:
    if ("transcript_support_level" %in% colnames(Transcripts)) {
        candidate.introns[Transcripts, on = "transcript_id",
            c("transcript_support_level") :=
                list(get("i.transcript_support_level"))
        ]
        candidate.introns[, c("transcript_support_level") :=
            tstrsplit(get("transcript_support_level"), split = " ")[[1]]]
        candidate.introns[
            is.na(get("transcript_support_level")),
            c("transcript_support_level") := "NA"
        ]
    }

    # Annotate protein_id:
    if ("protein_id" %in% colnames(Proteins)) {
        Proteins.red <- unique(Proteins[, c("transcript_id", "protein_id")])
        candidate.introns[Proteins.red,
            on = "transcript_id",
            c("protein_id") := list(get("i.protein_id"))
        ]
    }

    # Annotate ccds_id:
    if ("ccds_id" %in% colnames(Exons)) {
        Exons.red <- unique(Exons[, c("transcript_id", "ccds_id")])
        candidate.introns[Exons.red,
            on = "transcript_id",
            c("ccds_id") := list(get("i.ccds_id"))
        ]
    }
    return(candidate.introns)
}

#############################################################################
.process_introns_group <- function(candidate.introns,
        Exons_group.stranded, Exons_group.unstranded) {
    # Cannot annotate candidate introns by min and max exon_groups
    # because retained introns will overlap one or more exon groups
    # need to walk start -=1, end += 1, then do the overlap thing
    candidate.introns[, c("intron_start") := get("start")]
    candidate.introns[, c("intron_end") := get("end")]

    candidate.introns <- .process_introns_group_overlap(
        candidate.introns, Exons_group.stranded,
        c("gene_group_stranded", "exon_group_stranded_upstream",
            "gene_group_stranded", "exon_group_stranded_downstream"),
        c("gene_group", "exon_group", "gene_group", "exon_group")
    )
    # Need fix for retained_introns or sense_intronic where junction extends
    #   into the obligate introns
    candidate.introns <- .process_introns_group_fix_RI(
        candidate.introns, Exons_group.stranded,
        c("gene_group_stranded", "exon_group_stranded_upstream",
            "gene_group_stranded", "exon_group_stranded_downstream"),
        c("gene_group", "intron_number", "gene_group", "intron_number")
    )

    # Now repeat the same for unstranded condition
    candidate.introns <- .process_introns_group_overlap(
        candidate.introns, Exons_group.unstranded,
        c("gene_group_unstranded", "exon_group_unstranded_upstream",
            "gene_group_unstranded", "exon_group_unstranded_downstream"),
        c("gene_group", "exon_group", "gene_group", "exon_group")
    )

    # Need fix for retained_introns or sense_intronic where
    #   junction extends into the obligate introns
    candidate.introns <- .process_introns_group_fix_RI(
        candidate.introns, Exons_group.unstranded,
        c("gene_group_unstranded", "exon_group_unstranded_upstream",
            "gene_group_unstranded", "exon_group_unstranded_downstream"),
        c("gene_group", "intron_number", "gene_group", "intron_number")
    )

    # reset
    candidate.introns[, c("start") := get("intron_start")]
    candidate.introns[, c("end") := get("intron_end")]
    candidate.introns[, c("Event") := paste0(
        get("seqnames"), ":", get("intron_start"), "-",
        get("intron_end"), "/", get("strand")
    )]
    return(candidate.introns)
}

.process_introns_group_overlap <- function(target.DT, groups.DT,
    target.columns, groups.columns) {

    # Compile overlaps between upstream SS and exon groups definition
    OL <- .overlaps_exon_island(target.DT, groups.DT, upstream = TRUE)
    set(target.DT, from(OL), target.columns[1],
        groups.DT[, get(groups.columns[1])][to(OL)]
    ) # Gene group for upstream splice site
    set(target.DT, from(OL), target.columns[2],
        groups.DT[, get(groups.columns[2])][to(OL)]
    ) # Exon group for upstream splice site

    # Repeat for downstream SS
    OL <- .overlaps_exon_island(target.DT, groups.DT, upstream = FALSE)
    set(target.DT, from(OL),
        target.columns[3],
        groups.DT[, get(groups.columns[3])][to(OL)]
    )
    set(target.DT, from(OL),
        target.columns[4],
        groups.DT[, get(groups.columns[4])][to(OL)]
    )
    return(target.DT)
}

.process_introns_group_fix_RI <- function(
        target.DT, groups.DT,
        target.columns, groups.columns) {

    # Create meta-introns: introns formed between adjacent exon groups
    tmp <- .grDT(groups.DT, keep.extra.columns = TRUE)
    tmp.Introns_group <- .grlGaps(
        split(tmp, tmp$gene_group)
    )
    tmp.Introns_group <- as.data.table(tmp.Introns_group)
    setnames(tmp.Introns_group, "group_name", "gene_group")
    tmp.Introns_group[, c("intron_number") :=
        data.table::rowid(get("gene_group"))]
    tmp.Introns_group[get("strand") == "-",
        c("intron_number") :=
            max(get("intron_number")) + 1 - get("intron_number"),
        by = "gene_group"
    ]

    # Find introns that do not have annotated flanking exon islands
    # - these are typically introns of retained_intron transcripts as their
    #   junctions span into otherwise-obligate introns
    target.DT.subset <- target.DT[is.na(get(target.columns[2]))]
    target.DT <- target.DT[!is.na(get(target.columns[2]))]
    OL <- .overlaps_exon_island(
        target.DT.subset, tmp.Introns_group, upstream = TRUE
    )
    set(target.DT.subset, from(OL),
        target.columns[1],
        as.integer(tmp.Introns_group[, get(groups.columns[1])][to(OL)])
    ) # Gene group for upstream splice site
    set(target.DT.subset, from(OL),
        target.columns[2],
        tmp.Introns_group[, get(groups.columns[2])][to(OL)]
    ) # Exon group for upstream splice site
    target.DT <- rbind(target.DT, target.DT.subset)

    # Repeat for downstream SS
    target.DT.subset <- target.DT[is.na(get(target.columns[4]))]
    target.DT <- target.DT[!is.na(get(target.columns[4]))]
    OL <- .overlaps_exon_island(
        target.DT.subset,
        tmp.Introns_group,
        upstream = FALSE
    )
    set(target.DT.subset, from(OL),
        target.columns[3],
        as.integer(tmp.Introns_group[, get(groups.columns[3])][to(OL)])
    ) # Gene group for downstream splice site
    set(target.DT.subset, from(OL),
        target.columns[4],
        as.integer(tmp.Introns_group[, get(groups.columns[4])][to(OL)] + 1)
    ) # Exon group for downstream splice site
    return(rbind(target.DT, target.DT.subset))
}

# Creates GRanges that overlap upstream or downstream splice site by 1 nt
# - this allows assessment of which exon islands each intron bridges
#   by their overlap with flanking exon islands
.overlaps_exon_island <- function(intron.DT, groups.DT, upstream = TRUE) {
    if (all(c("intron_start", "intron_end") %in% colnames(intron.DT))) {
        int.DT <- intron.DT[, c("seqnames", "start", "end", "strand",
            "intron_start", "intron_end")]
    } else {
        int.DT <- intron.DT[, c("seqnames", "start", "end", "strand")]
        int.DT[, c("intron_start", "intron_end") :=
            list(get("start"), get("end"))]
    }
    if (upstream) {
        int.DT[get("strand") == "+", c("start") := get("intron_start") - 1]
        int.DT[get("strand") == "+", c("end") := get("intron_start")]
        int.DT[get("strand") == "-", c("start") := get("intron_end")]
        int.DT[get("strand") == "-", c("end") := get("intron_end") + 1]
    } else {
        int.DT[get("strand") == "-", c("start") := get("intron_start") - 1]
        int.DT[get("strand") == "-", c("end") := get("intron_start")]
        int.DT[get("strand") == "+", c("start") := get("intron_end")]
        int.DT[get("strand") == "+", c("end") := get("intron_end") + 1]
    }
    OL <- findOverlaps(
        .grDT(int.DT),
        .grDT(groups.DT)
    )
    return(OL)
}

################################################################################
# Sub

.gen_irf <- function(reference_path, extra_files, genome, chromosome_aliases) {
    .log("Generating SpliceWiz reference", "message")

    # Generating SpliceWiz-base references
    message("...prepping data")
    data <- .gen_irf_prep_data(reference_path)
    data2 <- .gen_irf_prep_introns(
        data[["candidate.introns"]], data[["Exons"]], extra_files)
    data2[["introns.unique"]] <- .gen_irf_prep_introns_unique(
        data2[["introns.unique"]], data2[["exclude.directional"]],
        data[["Genes.rev"]], data[["Genes.Extended"]]
    )
    message("...determining measurable introns (directional)")
    tmpdir.IntronCover.summa <- .gen_irf_export_introncover(
        .gen_irf_exclusion_zones(
            data2[["introns.unique"]], data2[["exclude.omnidirectional"]],
            data2[["exclude.directional"]], stranded = TRUE
        ), stranded = TRUE, reference_path, data2[["introns.unique"]]
    )
    message("...determining measurable introns (non-directional)")
    tmpnd.IntronCover.summa <- .gen_irf_export_introncover(
        tmpnd.IntronCover <- .gen_irf_exclusion_zones(
            data2[["introns.unique"]], data2[["exclude.omnidirectional"]],
            data2[["exclude.directional"]], stranded = FALSE
        ), stranded = FALSE, reference_path, data2[["introns.unique"]]
    )
    message("...writing ref-cover.bed")
    ref.cover <- .gen_irf_refcover(reference_path)
    message("...writing ref-ROI.bed")
    ref.ROI <- .gen_irf_ROI(reference_path, extra_files, genome,
        data[["Genes"]], data[["Transcripts"]])
    message("...writing ref-read-continues.ref")
    readcons <- .gen_irf_readcons(reference_path,
        tmpdir.IntronCover.summa, tmpnd.IntronCover.summa)
    message("...writing ref-sj.ref")
    ref.sj <- .gen_irf_sj(reference_path)

    chr <- data.frame(seqnames = names(GenomeInfoDb::seqinfo(genome)),
        seqlengths = unname(GenomeInfoDb::seqlengths(genome)))
    if (!is.null(chromosome_aliases)) {
        colnames(chromosome_aliases) <- c("name", "alias")
        chr$seqalias <- chromosome_aliases$alias[
            match(chr$seqnames, chromosome_aliases$name)]
        chr$seqalias[is.na(chr$seqalias)] <- ""
    } else {
        chr$seqalias <- ""
    }
    .gen_irf_final(reference_path, ref.cover, readcons, ref.ROI, ref.sj, chr)
    message("SpliceWiz reference generation completed")
}
################################################################################

# Load Genes, Exons, Transcripts, and Introns.
# Filter introns to protein_coding, processed_tx, lincRNA, antisense, and NMD
.gen_irf_prep_data <- function(reference_path) {
    Genes <- .grDT(
        read.fst(file.path(reference_path, "fst", "Genes.fst")),
        keep.extra.columns = TRUE
    )
    Genes.rev <- Genes
    strand(Genes.rev) <- ifelse(strand(Genes.rev) == "+", "-",
        ifelse(strand(Genes.rev) == "-", "+", "*")
    ) # Invert strand
    Genes.Extended <- reduce(c(
        flank(Genes.rev, 5000),
        flank(Genes.rev, 1000, start = FALSE)
    )) # 1000 nt upstream and 5000 nt downstream

    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )
    Exons <- .grDT(
        read.fst(file.path(reference_path, "fst", "Exons.fst")),
        keep.extra.columns = TRUE
    )
    Transcripts <- .grDT(
        read.fst(file.path(reference_path, "fst", "Transcripts.fst")),
        keep.extra.columns = TRUE
    )

    allowed_tx <- c("protein_coding", "processed_transcript",
            "lincRNA", "antisense", "nonsense_mediated_decay")
    candidate.introns <- candidate.introns[
        get("transcript_biotype") %in% allowed_tx]
    candidate.introns[, c("transcript_biotype") :=
        factor(get("transcript_biotype"), allowed_tx, ordered = TRUE)]

    # If common introns between several transcripts, sort by importance
    if ("transcript_support_level" %in% colnames(candidate.introns)) {
        # Sort by tsl first, then reverse later
        setorderv(candidate.introns, c("seqnames", "start", "end", "strand",
            "transcript_biotype", "transcript_support_level"))
    } else {
        setorderv(candidate.introns, c("seqnames", "start", "end", "strand",
            "transcript_biotype"))
    }

    final <- list(
        Genes = Genes, Genes.rev = Genes.rev,
        Genes.Extended = Genes.Extended,
        Exons = Exons, Transcripts = Transcripts,
        candidate.introns = candidate.introns
    )
    return(final)
}

.gen_irf_convert_seqnames <- function(gr, style) {
    if (!is(gr, "GRanges")) {
        .log("Cannot convert seqnames as object is not GRanges")
    }
    if (style != "") seqlevelsStyle(gr) <- style
    gr
}

# Unique introns, exclusion zones, known exon annotation
.gen_irf_prep_introns <- function(candidate.introns, Exons, extra_files) {

    # Filter for unique introns (same start / end)
    introns.unique <- unique(candidate.introns,
        by = c("seqnames", "start", "end", "width", "strand"))
    setorderv(introns.unique, c("seqnames", "start", "end", "strand"))
    introns.unique <- .grDT(introns.unique, keep.extra.columns = TRUE)

    # Directional exclusion - Regions overlapped by exons
    exclude.directional <- as.data.table(
        Exons[!grepl("intron", Exons$transcript_biotype)])
    exclude.directional <- unique(exclude.directional,
        by = c("seqnames", "start", "end", "width", "strand"))

    # Annotate known exons early.
    introns.unique.exon.dir <- findOverlaps(introns.unique,
        .grDT(exclude.directional), type = "within"
    )
    introns.unique.exon.nd <- findOverlaps(introns.unique,
        .grDT(exclude.directional), type = "within", ignore.strand = TRUE
    )

    # Mark introns if they exist within other exons (known RI events)
    introns.unique$known_exon_dir <-
        (seq_len(length(introns.unique)) %in% introns.unique.exon.dir@from)
    introns.unique$known_exon_nd <-
        (seq_len(length(introns.unique)) %in% introns.unique.exon.nd@from)

    # After known exon annotation, expand this to exclude 5 more nt's
    exclude.directional[, c("start", "end") :=
        list(get("start") - 5, get("end") + 5)]

    # Non-directional exclusion - Low mappability regions, blacklist regions
    exclude.omnidirectional <- GRanges(NULL)
    if (extra_files$MappabilityFile != "") {
        exclude.omnidirectional <- c(exclude.omnidirectional,
            .gen_irf_convert_seqnames(
                rtracklayer::import(extra_files$MappabilityFile, "bed"),
                extra_files$genome_style
            )
        )
    }
    if (extra_files$BlacklistFile != "") {
        exclude.omnidirectional <- c(exclude.omnidirectional,
            .gen_irf_convert_seqnames(
                rtracklayer::import(extra_files$BlacklistFile, "bed"),
                extra_files$genome_style
            )
        )
    }

    # merge with any gaps <= 9
    exclude.omnidirectional <-
        reduce(exclude.omnidirectional, min.gapwidth = 9)

    # Filter out introns are lie completely within low mappability or blacklists
    if (length(exclude.omnidirectional) > 0) {
        introns.unique.blacklisted <- findOverlaps(introns.unique,
            exclude.omnidirectional, type = "within"
        )
        introns.unique <- introns.unique[-introns.unique.blacklisted@from]
    }

    final <- list(
        introns.unique = introns.unique,
        exclude.directional = exclude.directional,
        exclude.omnidirectional = exclude.omnidirectional
    )
    return(final)
}

# Annotate anti-near, anti-over
.gen_irf_prep_introns_unique <- function(introns.unique, exclude.directional,
        Genes.rev, Genes.Extended) {

    # Antiover: overlaps within anti-sense genes
    # Antinear: overlaps within 1000 / 5000 nt up/downstream of antisense gene
    introns.unique.antiover <- findOverlaps(introns.unique, Genes.rev)
    introns.unique.antinear <- findOverlaps(introns.unique, Genes.Extended)

    introns.unique$antiover <-
        (seq_len(length(introns.unique)) %in% introns.unique.antiover@from)
    introns.unique$antinear <-
        (seq_len(length(introns.unique)) %in% introns.unique.antinear@from)

    introns.unique$intron_width <- BiocGenerics::width(introns.unique)

    # Remove introns less than 50 bp:
    introns.unique <- introns.unique[BiocGenerics::width(introns.unique) > 50]

    # remove 5 bases from start & end
    BiocGenerics::start(introns.unique) <-
        BiocGenerics::start(introns.unique) + 5
    BiocGenerics::end(introns.unique) <-
        BiocGenerics::end(introns.unique) - 5

    # NB intron_start and intron_end represent actual start and end of introns
    return(introns.unique)
}

#
.gen_irf_exclusion_zones <- function(introns.unique,
        exclude.omnidirectional, exclude.directional,
        stranded = TRUE) {
    if (!stranded) {
        exclude.directional.reverse <- copy(exclude.directional)
        exclude.directional.reverse[get("strand") == "-", c("strand") := "P"]
        exclude.directional.reverse[get("strand") == "+", c("strand") := "-"]
        exclude.directional.reverse[get("strand") == "P", c("strand") := "+"]
        exclude.directional.gr <- c(.grDT(exclude.directional),
            .grDT(exclude.directional.reverse))
    } else {
        exclude.directional.gr <- .grDT(exclude.directional)
    }
    if (length(exclude.omnidirectional) > 0) {
        introns.intersect <- GenomicRanges::intersect(introns.unique,
            c(exclude.omnidirectional, exclude.directional.gr))
    } else if (length(exclude.directional) > 0) {
        introns.intersect <- GenomicRanges::intersect(introns.unique,
            exclude.directional.gr)
    }
    # introns.intersect is the list of intron regions that
    #   should be excluded from analysis

    OL <- findOverlaps(introns.unique, introns.intersect)
    # make a GRanges same size as the number of intersections
    introns.intersect.final <- introns.intersect[to(OL)]
    introns.intersect.final$intron_id <- introns.unique$intron_id[from(OL)]

    # Create GRangesList of intron regions and exclusions, split by intron_id
    introns.unique.ID <- split(introns.unique, introns.unique$intron_id)
    introns.intersect.ID <- split(introns.intersect.final,
            introns.intersect.final$intron_id)
    # Filter for introns that have exclusion regions
    introns.unique.ID.compare <- introns.unique.ID[
        names(introns.unique.ID) %in% names(introns.intersect.ID)
    ]

    # Take the difference: returns intron regions that need to be measured
    #   GRangesList split by intron_id
    IntronCover <- setdiff(introns.unique.ID.compare, introns.intersect.ID)

    # Now add back introns that did not require intersection
    #   (or would have been excluded as known-exons)
    IntronCover <- c(IntronCover,
        introns.unique.ID[!(names(introns.unique.ID) %in%
            names(introns.intersect.ID))])
    if (stranded) {
        IntronCover <- c(IntronCover,
            introns.unique.ID[names(introns.unique.ID) %in%
                introns.unique$intron_id[introns.unique$known_exon_dir == TRUE]]
        )
    } else {
        IntronCover <- c(IntronCover,
            introns.unique.ID[names(introns.unique.ID) %in%
                introns.unique$intron_id[introns.unique$known_exon_nd == TRUE]]
        )
    }
    IntronCover <- as.data.table(IntronCover)
    IntronCover <- IntronCover[,
        c("seqnames", "start", "end", "strand", "width", "group_name")
    ]
    colnames(IntronCover)[6] <- "intron_id"
    return(IntronCover)
}

# Annotates IntronCover data
.gen_irf_export_introncover <- function(IntronCover, stranded = TRUE,
        reference_path, introns.unique) {
    IntronCover.summa <- IntronCover
    IntronCover.summa[,
        c("num_blocks", "inclbases") := list(.N, sum(get("width"))),
        by = "intron_id"]
    IntronCover.summa <- unique(
        IntronCover.summa[, c("intron_id", "num_blocks", "inclbases")],
        by = "intron_id")
    IntronCover.summa[as.data.table(introns.unique),
        on = "intron_id",
        c("seqnames", "intron_start", "intron_end",
            "intron_width", "width", "strand", "gene_name", "transcript_id")
        := list(get("i.seqnames"), get("i.intron_start"),
                get("i.intron_end"), get("i.intron_width"),
                get("i.width"), get("i.strand"), get("i.gene_name"),
                get("i.transcript_id"))
    ]
    if (stranded) {
        IntronCover.summa[as.data.table(introns.unique), on = "intron_id",
            c("known_exon_dir", "GG", "EG_up", "EG_down") :=
            list(get("i.known_exon_dir"), get("i.gene_group_stranded"),
                get("i.exon_group_stranded_upstream"),
                get("i.exon_group_stranded_downstream"))
        ]
    } else {
        IntronCover.summa[as.data.table(introns.unique), on = "intron_id",
            c("known_exon_nd", "antiover", "antinear",
                "GG", "EG_up", "EG_down") :=
            list(get("i.known_exon_nd"), get("i.antiover"), get("i.antinear"),
                get("i.gene_group_unstranded"),
                get("i.exon_group_unstranded_upstream"),
                get("i.exon_group_unstranded_downstream"))
        ]
    }
    IntronCover.summa[, c("exclbases") :=
        get("intron_width") - get("inclbases")]
    # Exclude exclbases / width > 0.3
    IntronCover.summa <-
        IntronCover.summa[get("exclbases") / get("intron_width") < 0.3]
    IntronCover <- .semi_join_DT(IntronCover, IntronCover.summa,
        by = "intron_id")
    IntronCover.summa <- .gen_irf_irfname(IntronCover.summa,
        stranded = stranded)

    IntronCover <- .grDT(IntronCover, keep.extra.columns = TRUE)
    IntronCover <- split(IntronCover, IntronCover$intron_id)
    names(IntronCover) <- IntronCover.summa$IRname[match(
        names(IntronCover), IntronCover.summa$intron_id)]

    # Arrange by seqnames, start, end, strand
    setorderv(IntronCover.summa,
        c("seqnames", "intron_start", "intron_end", "strand"))
    IntronCover <- IntronCover[IntronCover.summa$IRname]

    # Export as 12-column BED file
    rtracklayer::export(IntronCover, file.path(reference_path,
        ifelse(stranded, "tmpdir.IntronCover.bed", "tmpnd.IntronCover.bed")
    ))
    write.fst(IntronCover.summa, file.path(
        reference_path, "fst",
        ifelse(stranded, "Introns.Dir.fst", "Introns.ND.fst")
    ))
    return(IntronCover.summa)
}

# Generates SpliceWiz intron name
.gen_irf_irfname <- function(IntronCover.summa, stranded = TRUE) {
    if (stranded) {
        IntronCover.summa[, c("IRname") := paste("dir",
            get("gene_name"), get("intron_id"), get("strand"),
            get("num_blocks"), sprintf("%.f", get("intron_start") - 1),
            sprintf("%.f", get("intron_end")), get("inclbases"),
            get("exclbases"),
            ifelse(get("known_exon_dir"), "known-exon", "clean"), sep = "/"
        )]
    } else {
        IntronCover.summa[, c("IRname") := paste("nd",
            get("gene_name"), get("intron_id"), get("strand"),
            get("num_blocks"), sprintf("%.f", get("intron_start") - 1),
            sprintf("%.f", get("intron_end")), get("inclbases"),
            get("exclbases"), sep = "/"
        )]
        # casewise naming of last condition
        IntronCover.summa[
            get("known_exon_nd") & get("antiover") & get("antinear"),
            c("IRname") := paste(get("IRname"),
                "known-exon+anti-over+anti-near", sep = "/")]
        IntronCover.summa[
            get("known_exon_nd") & get("antiover") & !get("antinear"),
            c("IRname") := paste(get("IRname"),
                "known-exon+anti-over", sep = "/")]
        IntronCover.summa[
            get("known_exon_nd") & !get("antiover") & get("antinear"),
            c("IRname") := paste(get("IRname"),
                "known-exon+anti-near", sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & get("antiover") & get("antinear"),
            c("IRname") := paste(get("IRname"),
                "anti-over+anti-near", sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & !get("antiover") & get("antinear"),
            c("IRname") := paste(get("IRname"), "anti-near", sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & get("antiover") & !get("antinear"),
            c("IRname") := paste(get("IRname"), "anti-over", sep = "/")]
        IntronCover.summa[
            get("known_exon_nd") & !get("antiover") & !get("antinear"),
            c("IRname") := paste(get("IRname"), "known-exon", sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & !get("antiover") & !get("antinear"),
            c("IRname") := paste(get("IRname"), "clean", sep = "/")]
    }
    return(IntronCover.summa)
}

# Imports IntronCover temp files, generates SpliceWiz-readable reference
.gen_irf_refcover <- function(reference_path) {
    tmpdir.IntronCover <- fread(file.path(
        reference_path, "tmpdir.IntronCover.bed"
    ), sep = "\t")
    tmpdir.IntronCover[, c("cat") := "dir"]
    tmpnd.IntronCover <- fread(file.path(
        reference_path, "tmpnd.IntronCover.bed"
    ), sep = "\t")
    tmpnd.IntronCover[, c("cat") := "nd"]

    ref.cover <- rbind(tmpdir.IntronCover, tmpnd.IntronCover)
    setorderv(ref.cover, c("V1", "V2", "V3", "V6", "cat"))
    ref.cover$cat <- NULL
    ref.cover[, c("V9") := as.character(get("V9"))]
    ref.cover[, c("V9") := "255,0,0"]

    return(ref.cover)
}

.gen_irf_ROI <- function(reference_path, extra_files, genome,
        Genes, Transcripts) {

    # List of rRNA regions
    rRNA <- as.data.table(Transcripts[grepl("rRNA", Transcripts$gene_biotype)])
    if (nrow(rRNA) > 0) {
        rRNA[, c("start") := get("start") - 1]
        rRNA[, c("name") := paste("rRNA", get("seqnames"), get("start"),
            get("end"), get("strand"), get("transcript_id"),
            get("gene_biotype"), get("gene_id"), get("gene_name"), sep = "/")]
        rRNA <- rRNA[, c("seqnames", "start", "end", "name"), with = FALSE]
    } else {
        rRNA <- c()
    }

    # List of nonPolyA regions
    if (extra_files$nonPolyAFile != "") {
        nonPolyA <- .gen_irf_convert_seqnames(
            rtracklayer::import(extra_files$nonPolyAFile, "bed"),
            extra_files$genome_style
        )

        nonPolyA <- as.data.table(nonPolyA)
        nonPolyA <- nonPolyA[, c("seqnames", "start", "end"), with = FALSE]
        nonPolyA[, c("name") := "NonPolyA"]
    } else {
        nonPolyA <- c()
    }

    # Intergenic regions
    AllChr <- makeGRangesListFromDataFrame(data.frame(
        seqnames = names(seqinfo(genome)),
        start = 1, end = seqlengths(seqinfo(genome)),
        names = names(seqinfo(genome))
    ), split.field = "names")
    Genes.chr <- c(
        Genes, flank(Genes, 10000), flank(Genes, 10000, start = FALSE)
    )
    Genes.chr <- reduce(Genes.chr, min.gapwidth = 1000)
    Genes.chr$chr <- seqnames(Genes.chr)
    Genes.chr <- split(Genes.chr, Genes.chr$chr)
    AllChr <- AllChr[names(Genes.chr)]
    AllChr.split <- setdiff(AllChr, Genes.chr, ignore.strand = TRUE)
    Intergenic <- unlist(AllChr.split)
    if (length(Intergenic) > 0) {
        names(Intergenic) <- seq_len(length(Intergenic))
        Intergenic <- as.data.table(Intergenic)
        Intergenic <- Intergenic[, c("seqnames", "start", "end"), with = FALSE]
        Intergenic[, c("name") := paste("Intergenic", Intergenic$seqnames,
            sep = "/")]
    } else {
        Intergenic <- c()
    }
    ref.ROI <- rbind(rRNA, nonPolyA, Intergenic)
    if (!is.null(ref.ROI) && nrow(ref.ROI) > 0) {
        # ref.ROI = as.data.table(ref.ROI)
        setorderv(ref.ROI, c("seqnames", "start"))
        ref.ROI[, c("start") := get("start") - 1] # convert back to 0-based
    }
    return(ref.ROI)
}

.gen_irf_readcons <- function(reference_path,
        tmpdir.IntronCover.summa, tmpnd.IntronCover.summa) {
    # ref-read-continues.ref
    introns.unique.readcons <- rbind(
        tmpdir.IntronCover.summa[,
            c("seqnames", "intron_start", "intron_end", "strand")],
        tmpnd.IntronCover.summa[,
            c("seqnames", "intron_start", "intron_end", "strand")]
    )
    # 0-based
    introns.unique.readcons[, c("intron_start") := get("intron_start") - 1]

    readcons.left <- introns.unique.readcons[
        ,
        c("seqnames", "intron_start", "strand")
    ]
    readcons.right <- introns.unique.readcons[
        ,
        c("seqnames", "intron_end", "strand")
    ]
    colnames(readcons.left) <- c("V1", "V2", "V3")
    colnames(readcons.right) <- c("V1", "V2", "V3")
    readcons <- rbind(readcons.left, readcons.right)
    setorderv(readcons, c("V1", "V2", "V3"))
    readcons <- unique(readcons)
    gc()
    return(readcons)
}

.gen_irf_sj <- function(reference_path) {
    # ref-sj.ref
    # Reload candidate introns here, as we've filtered this before
    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )

    ref.sj <- candidate.introns[, c("seqnames", "start", "end", "strand")]
    ref.sj$Is_NMD <- ifelse(
        grepl("nonsense_mediated_decay", candidate.introns$transcript_biotype),
        "NMD", ""
    )

    # annotate NMD-unique junctions
    ref.sj <- ref.sj[, lapply(.SD, function(x) {
        ifelse(all(x != ""), "NMD", "")
    }), by = c("seqnames", "start", "end", "strand")]

    # BED file conversion
    ref.sj[, c("start") := get("start") - 1]
    setorderv(ref.sj, c("seqnames", "start", "end", "strand"))
    gc()
    return(ref.sj)
}

.gen_irf_final <- function(reference_path,
        ref.cover, readcons, ref.ROI, ref.sj,
        chromosome_aliases
) {
    IRF_file <- file.path(reference_path, "SpliceWiz.ref")
    # Concatenate all 4 reference files into one file
    fwrite(list("# ref-cover.bed"), IRF_file,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    fwrite(ref.cover, IRF_file, append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    fwrite(list("# ref-read-continues.ref"), IRF_file, append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    fwrite(readcons, IRF_file, append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    fwrite(list("# ref-ROI.bed"), IRF_file, append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    if (!is.null(ref.ROI) && nrow(ref.ROI) > 0) {
        fwrite(ref.ROI, IRF_file, append = TRUE,
            sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    }
    fwrite(list("# ref-sj.ref"), IRF_file, append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    fwrite(ref.sj, IRF_file, append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)

    if (!is.null(chromosome_aliases)) {
        fwrite(list("# ref-chrs.ref"), IRF_file, append = TRUE,
            sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
        fwrite(chromosome_aliases, IRF_file, append = TRUE,
            sep = "\t", eol = "\n", col.names = FALSE, scipen = 50)
    }

    # Add EOF (to avoid undefined behaviour when there is no termination char)
    fwrite(list("# EOF"), IRF_file, append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )

    gzip(filename = IRF_file, destname = paste0(IRF_file, ".gz"),
        overwrite = TRUE)
    if (file.exists(IRF_file) & file.exists(paste0(IRF_file, ".gz"))) {
        file.remove(IRF_file)
    }
    # cleanup
    if (file.exists(file.path(reference_path, "tmpdir.IntronCover.bed"))) {
        file.remove(file.path(reference_path, "tmpdir.IntronCover.bed"))
    }
    if (file.exists(file.path(reference_path, "tmpnd.IntronCover.bed"))) {
        file.remove(file.path(reference_path, "tmpnd.IntronCover.bed"))
    }
}
################################################################################

# Determines which spliced / IR transcripts are NMD substrates
# Assumes NMD substrates if PTC is < 50 nt from last EJC
.gen_nmd <- function(reference_path, genome) {

    Exons.tr <- .gen_nmd_exons_trimmed(reference_path)
    protein.introns <- .gen_nmd_protein_introns(reference_path, Exons.tr)

    NMD.Table <- .gen_nmd_determine(Exons.tr, protein.introns, genome, 50)
    protein.introns.red <- unique(
        protein.introns[, c("intron_id", "intron_type")])
    NMD.Table[protein.introns.red, on = "intron_id",
        c("intron_type") := get("i.intron_type")]

    write.fst(NMD.Table, file.path(reference_path, "fst", "IR.NMD.fst"))
    gc()
}

# Get exons but trim by start codon
.gen_nmd_exons_trimmed <- function(reference_path) {
    Exons.tr <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.fst"))
    )
    Misc <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Misc.fst"))
    )
    start.DT <- Misc[get("type") == "start_codon"]
    Exons.tr <- Exons.tr[get("transcript_id") %in%
        start.DT[, get("transcript_id")]]

    Exons.tr[start.DT,
        on = c("transcript_id"),
        c("sc_start", "sc_end") := list(get("i.start"), get("i.end"))
    ]
    Exons.tr[
        get("start") < get("sc_start") & get("strand") == "+",
        c("start") := get("sc_start")
    ]
    Exons.tr[
        get("end") < get("sc_start") & get("strand") == "+",
        c("end") := get("sc_start")
    ]
    Exons.tr[
        get("start") > get("sc_end") & get("strand") == "-",
        c("start") := get("sc_end")
    ]
    Exons.tr[
        get("end") > get("sc_end") & get("strand") == "-",
        c("end") := get("sc_end")
    ]
    Exons.tr <- Exons.tr[get("start") < get("end")]
    return(Exons.tr)
}

# Get introns and annotate by whether they are 5', 3' or CDS
.gen_nmd_protein_introns <- function(reference_path, Exons.tr) {
    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )
    Misc <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Misc.fst"))
    )
    protein.introns <- candidate.introns[
        get("transcript_id") %in% Exons.tr$transcript_id
    ]
    # determine here whether protein introns are CDS, 5' or 3' UTR introns
    UTR5 <- Misc[get("type") == "five_prime_utr"]
    UTR5.introns <- .grlGaps(split(
        makeGRangesFromDataFrame(as.data.frame(UTR5)),
        UTR5$transcript_id
    ))
    UTR5.introns <- as.data.table(UTR5.introns)
    UTR3 <- Misc[get("type") == "three_prime_utr"]
    UTR3.introns <- .grlGaps(split(
        makeGRangesFromDataFrame(as.data.frame(UTR3)),
        UTR3$transcript_id
    ))
    UTR3.introns <- as.data.table(UTR3.introns)

    CDS.introns <- .grlGaps(split(
        makeGRangesFromDataFrame(as.data.frame(Exons.tr)),
        Exons.tr$transcript_id
    ))
    CDS.introns <- as.data.table(CDS.introns)

    protein.introns[UTR5.introns,
        on = c("seqnames", "start", "end", "strand"),
        c("intron_type") := "UTR5"
    ]
    protein.introns[UTR3.introns,
        on = c("seqnames", "start", "end", "strand"),
        c("intron_type") := "UTR3"
    ]
    protein.introns[CDS.introns,
        on = c("seqnames", "start", "end", "strand"),
        c("intron_type") := "CDS"
    ]
    return(protein.introns)
}

# Given a list of exons and introns, and genome sequence
# - Generate a list of whether spliced or unspliced transcripts are NMD subs
.gen_nmd_determine <- function(exon.DT, intron.DT, genome, threshold = 50) {
    .log("Predicting NMD transcripts from genome sequence", "message")
    exon.DT <- exon.DT[,
        c("seqnames", "start", "end", "strand", "transcript_id")]
    exon_gr <- .grDT(exon.DT)
    message("...exonic transcripts")

    set(exon.DT, , "seq", as.character(getSeq(genome, exon_gr)))
    final <- .gen_nmd_determine_spliced_exon(exon.DT, intron.DT,
        threshold = threshold)

    intron.DT.use <- intron.DT[get("intron_type") != "UTR5"]
    exon.DT.skinny <- exon.DT[, -("seq")]
    i_partition <- c(seq(1, nrow(intron.DT.use), by = 10000),
        nrow(intron.DT.use) + 1)
    message("...retained introns")
    pb <- txtProgressBar(max = length(i_partition) - 1, style = 3)
    l_seq <- 1000
    for (i in seq_len(length(i_partition) - 1)) {
        setTxtProgressBar(pb, i)
        dash_progress("Determining NMD Transcripts: Calculating...",
            length(i_partition) - 1)
        intron.part <- intron.DT.use[
            seq(i_partition[i], i_partition[i + 1] - 1),
            c("transcript_id", "intron_id",
                "seqnames", "start", "end", "strand")
        ]
        set(intron.part, , "type", "intron")

        intron.part.upstream <- .gen_nmd_determine_build_introns_upstream(
            intron.part, exon.DT.skinny, use_short = TRUE)
        intron.part.upstream <- .gen_nmd_determine_retrieve_short_seq(
            exon.DT, intron.part.upstream, genome, l_seq = l_seq)
        final <- .gen_nmd_determine_translate(
            final, intron.part.upstream, use_short = TRUE,
            threshold = threshold)
        intron_id_exclude <- unique(intron.part.upstream$intron_id)
        intron_id_exclude <- intron_id_exclude[(
            intron_id_exclude %in%
            final[get("IRT_is_NMD") == TRUE, get("intron_id")]
        )]
        intron.part <- intron.part[!(get("intron_id") %in% intron_id_exclude)]

        intron.part.upstream <- .gen_nmd_determine_build_introns_upstream(
            intron.part, exon.DT.skinny, use_short = FALSE)
        intron.part.upstream <- .gen_nmd_determine_retrieve_full_seq(
            exon.DT, intron.part.upstream, genome)
        final <- .gen_nmd_determine_translate(
            final, intron.part.upstream, use_short = FALSE,
            threshold = threshold)
    }
    setTxtProgressBar(pb, i)
    close(pb)
    message("done")
    return(final)
}

# Determine PTC pos, and NMD status of spliced transcripts
.gen_nmd_determine_spliced_exon <- function(
        exon.DT, intron.DT, threshold = 50) {
    exon.MLE.DT <- copy(exon.DT)
    setorderv(exon.MLE.DT, "start")
    exon.MLE.DT[, c("elem_number") := data.table::rowid(get("transcript_id"))]
    exon.MLE.DT[get("strand") == "-",
        c("elem_number") := max(get("elem_number")) + 1 - get("elem_number"),
        by = "transcript_id"]
    exon.MLE.DT[, by = "transcript_id",
        c("is_last_elem") := (get("elem_number") == max(get("elem_number")))]
    exon.MLE.DT <- exon.MLE.DT[get("is_last_elem") == FALSE]

    # sort by order
    setorderv(exon.MLE.DT, c("transcript_id", "elem_number"))
    exon.MLE.DT <- exon.MLE.DT[, c("transcript_id", "seq")]
    splice <- exon.MLE.DT[, lapply(.SD, paste0, collapse = ""),
        by = "transcript_id"]
    splice[as.numeric(regexpr("N", get("seq"))) < 0,
        c("AA") := as.character(
            Biostrings::translate(as(.trim_3(get("seq")), "DNAStringSet"))
        )
    ]
    # Find nucleotide position of first stop codon
    splice[, c("stop_pos") :=
        as.numeric(regexpr("\\*", get("AA"))) * 3 - 2]
    splice[get("stop_pos") < 0, c("stop_pos") := NA]
    splice[, c("splice_len") := nchar(get("seq"))]
    splice[!is.na(get("AA")),
        c("stop_to_EJ") := get("splice_len") - get("stop_pos")]
    intron.DT <- intron.DT[,
        c("seqnames", "start", "end", "strand", "transcript_id", "intron_id")]
    final <- intron.DT[, c("intron_id", "transcript_id")]
    final[splice[, c("transcript_id", "stop_pos", "splice_len", "stop_to_EJ")],
        on = "transcript_id",
        c("splice_stop_pos", "splice_start_to_last_EJ",
                "splice_stop_to_last_EJ") :=
            list(get("i.stop_pos"), get("i.splice_len"), get("i.stop_to_EJ"))
    ]
    final[, c("splice_is_NMD") :=
        ifelse(get("splice_start_to_last_EJ") - get("splice_stop_to_last_EJ")
        >= threshold, TRUE, FALSE)]
    final[is.na(get("splice_stop_to_last_EJ")), c("splice_is_NMD") := FALSE]
    final[is.na(get("splice_start_to_last_EJ")), c("splice_is_NMD") := NA]
    return(final)
}

# Builds transcripts with one retained intron
# Then trims the last exon so the transcript terminus is the last EJC
.gen_nmd_determine_build_introns_upstream <- function(
        intron.part, exon.DT.skinny, use_short = FALSE
) {
    intron.part.skinny <- intron.part[, c("transcript_id", "intron_id")]
    # join exons with introns to determine phase of intron
    exon.DT.skinny.copy <- copy(exon.DT.skinny)
    exon.DT.skinny.copy <- exon.DT.skinny.copy[
        intron.part.skinny, on = "transcript_id",
        allow.cartesian = TRUE]
    exon.DT.skinny.copy <- exon.DT.skinny.copy[,
        c("transcript_id", "intron_id", "seqnames", "start", "end", "strand")]
    set(exon.DT.skinny.copy, , "type", "exon")
    intron.part.upstream <- rbindlist(list(intron.part, exon.DT.skinny.copy))

    setorderv(intron.part.upstream, c("seqnames", "start"))
    intron.part.upstream[,
        c("elem_number") := data.table::rowid(get("intron_id"))]
    intron.part.upstream[get("strand") == "-",
        c("elem_number") :=
            max(get("elem_number")) + 1 - get("elem_number"),
        by = "intron_id"
    ]

    # trim exons downstream of intron
    intron.part.upstream.intron <- intron.part.upstream[get("type") == "intron"]
    intron.part.upstream[intron.part.upstream.intron,
        on = "intron_id", c("intron_pos") := get("i.elem_number")
    ]
    intron.part.upstream <- intron.part.upstream[!is.na(get("intron_pos"))]
    if (use_short) {
        intron.part.upstream <-
            intron.part.upstream[get("elem_number") < get("intron_pos") |
                get("type") == "intron"]
    } else {
        # remove last exon: then the terminus is the last exon junction
        intron.part.upstream[,
            by = "transcript_id",
            c("is_last_elem") := (get("elem_number") == max(get("elem_number")))
        ]
        intron.part.upstream <-
            intron.part.upstream[!get("is_last_elem") |
                get("type") == "intron"]

    }
    return(intron.part.upstream)
}

# Retrieves full sequence of IR-transcript, up to last EJC
.gen_nmd_determine_retrieve_full_seq <- function(
        exon.DT, intron.part.upstream, genome
) {
    intron.part.short <- intron.part.upstream[get("type") == "intron"]
    intron.short_gr <- .grDT(intron.part.short)
    intron.part.short[,
        c("seq") := as.character(getSeq(genome, intron.short_gr))]
    intron.part.upstream[exon.DT,
        on = c("transcript_id", "seqnames", "start", "end", "strand"),
        c("seq") := get("i.seq")
    ]
    intron.part.upstream[intron.part.short,
        on = c("intron_id", "type"),
        c("seq") := get("i.seq")
    ]
    return(intron.part.upstream)
}

# Retrieve sequence up to `l_seq` bases of retained intron
# Saves having to generate full sequence if the PTC is early
.gen_nmd_determine_retrieve_short_seq <- function(
        exon.DT, intron.part.upstream, genome, l_seq = 1000, threshold = 50
) {
    intron.part.short <- intron.part.upstream[get("type") == "intron"]
    # Truncate intron by threshold as its terminus is taken as (EJC - threshold)
    intron.part.short[get("strand") == "+" &
        get("end") - get("start") > threshold,
        c("end") := get("end") - threshold]
    intron.part.short[get("strand") == "-" &
        get("end") - get("start") > threshold,
        c("start") := get("start") + threshold]
    intron.part.short[
        get("strand") == "+" & get("end") - get("start") > l_seq,
        c("end") := get("start") + l_seq
    ]
    intron.part.short[
        get("strand") == "-" & get("end") - get("start") > l_seq,
        c("start") := get("end") - l_seq
    ]

    intron.short_gr <- .grDT(intron.part.short)
    intron.part.short[,
        c("seq") := as.character(getSeq(genome, intron.short_gr))]
    intron.part.upstream[exon.DT,
        on = c("transcript_id", "seqnames", "start", "end", "strand"),
        c("seq") := get("i.seq")
    ]
    intron.part.upstream[intron.part.short,
        on = c("intron_id", "type"),
        c("seq") := get("i.seq")
    ]
    return(intron.part.upstream)
}

# Translate nucleotide sequence to determine PTC position
.gen_nmd_determine_translate <- function(
        splice_table, elems, use_short = FALSE, threshold = 50
) {
    setorderv(elems, c("transcript_id", "elem_number"))
    elems <- elems[, c("intron_id", "seq")]

    IRT <- elems[,
        lapply(.SD, paste0, collapse = ""), by = "intron_id"]
    # trim
    IRT[, c("seq") := substr(get("seq"), 1,
        nchar(get("seq")) - (nchar(get("seq")) %% 3))]
    IRT[
        as.numeric(regexpr("N", get("seq"))) < 0,
        c("AA") := as.character(
            Biostrings::translate(as(.trim_3(get("seq")), "DNAStringSet"))
        )
    ]

    # Find nucleotide position of first stop codon
    IRT[, c("stop_pos") :=
        as.numeric(regexpr("\\*", get("AA"))) * 3 - 2]
    IRT[get("stop_pos") < 0, c("stop_pos") := NA]
    IRT[, c("IRT_len") := nchar(get("seq"))]
    IRT[!is.na(get("AA")), c("stop_to_EJ") :=
        get("IRT_len") - get("stop_pos")]
    IRT[, c("use_short") := use_short]

    IRT[, c("IRT_is_NMD") := ifelse(
        get("stop_to_EJ") >= get("threshold"), TRUE, FALSE)]
    IRT[is.na(get("stop_pos")), c("IRT_is_NMD") := FALSE]
    IRT[is.na(get("IRT_len")), c("IRT_is_NMD") := NA]

    # Annotate into splice_table
    splice_table[IRT,
        on = "intron_id",
        c(
            "IRT_stop_pos", "IRT_start_to_last_EJ", "IRT_stop_to_last_EJ",
            "IRT_use_short", "IRT_is_NMD"
        ) :=
            list(
                get("i.stop_pos"), get("i.IRT_len"), get("i.stop_to_EJ"),
                get("i.use_short"), get("i.IRT_is_NMD")
            )
    ]
    return(splice_table)
}

################################################################################
# Sub

# Generate a list of ASEs
.gen_splice <- function(reference_path) {
    .log("Annotating Splice Events", "message")
    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )
    introns.skipcoord <- .gen_splice_skipcoord(
        reference_path, candidate.introns)

    message("Annotating Mutually-Exclusive-Exon Splice Events...",
        appendLF = FALSE
    )
    introns_found_MXE <- .gen_splice_MXE(introns.skipcoord)
    message("done")

    # annotate skipped junctions with two included junctions
    message("Annotating Skipped-Exon Splice Events...", appendLF = FALSE)
    introns_found_SE <- .gen_splice_SE(introns.skipcoord, candidate.introns)
    message("done")

    message("Annotating Alternate 5' / 3' Splice Site Splice Events...",
        appendLF = FALSE)

    introns_found_A5SS <- .gen_splice_A5SS(candidate.introns)
    introns_found_A3SS <- .gen_splice_A3SS(candidate.introns)
    message("done")

    message("Annotating Alternate First / Last Exon Splice Events...",
        appendLF = FALSE)
    # AFE/ALE

    introns_found_AFE <- .gen_splice_AFE(candidate.introns, introns_found_A5SS)
    introns_found_ALE <- .gen_splice_ALE(candidate.introns, introns_found_A3SS)
    message("done")

    # Annotate known RI's
    message("Annotating known retained introns...",
        appendLF = FALSE)
    introns_found_RI <- .gen_splice_RI(candidate.introns, reference_path)
    message("done")
    gc()

    #   Filter for valid splicing
    is_valid_splice_type <- function(x) !is.null(x) && nrow(x) > 0
    tmp_AS <- list(
        introns_found_MXE, introns_found_SE,
        introns_found_AFE, introns_found_ALE,
        introns_found_A5SS, introns_found_A3SS,
        introns_found_RI
    )
    tmp_AS <- base::Filter(is_valid_splice_type, tmp_AS)
    AS_Table <- rbindlist(tmp_AS)

    if (nrow(AS_Table) > 0) {
        .gen_splice_save(AS_Table, candidate.introns, reference_path)
        .log("Splice Annotations Filtered", "message")
    } else {
        message("No splice events found\n")
    }
}

################################################################################

# Generate a list of skip-coordinates
# - These are introns if the downstream exon is skipped
.gen_splice_skipcoord <- function(reference_path, candidate.introns) {
    GeneOrder <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Genes.fst"))
    )
    setorder(GeneOrder, seqnames, start, end, strand)

    introns.skipcoord <- copy(candidate.introns)
    introns.skipcoord[, c("gene_id") :=
        factor(get("gene_id"), GeneOrder[, get("gene_id")], ordered = TRUE)]

    setorderv(introns.skipcoord,
        c("gene_id", "transcript_name", "intron_number"))

    introns.skipcoord[, c("skip_coord") := ""]
    introns.skipcoord[get("strand") == "+",
        c("skip_coord") := paste0(
                get("seqnames"), ":", get("intron_start"),
                "-", data.table::shift(get("intron_end"), 1, NA, "lead"),
                "/", get("strand")
        ),
        by = "transcript_id"
    ]
    introns.skipcoord[get("strand") == "-",
        c("skip_coord") := paste0(
                get("seqnames"), ":",
                data.table::shift(get("intron_start"), 1, NA, "lead"),
                "-", get("intron_end"), "/", get("strand")
        ),
        by = "transcript_id"
    ]
    introns.skipcoord[grepl("NA", get("skip_coord")), c("skip_coord") := NA]

    introns.skipcoord[, c("skip_coord_2") :=
        data.table::shift(get("skip_coord"), 1, NA, "lag")]
    return(introns.skipcoord)
}

# Generate a list of MXE
.gen_splice_MXE <- function(introns.skipcoord) {
    introns_search_MXE <- introns.skipcoord[introns.skipcoord[,
        .I[get("intron_number") < max(get("intron_number"))],
        by = "transcript_id"]$V1]
    introns_search_MXE <- introns_search_MXE[
        introns_search_MXE[, .N, by = c("skip_coord")],
        on = c("skip_coord"),
        c("N") := get("i.N")]
    introns_search_MXE <- introns_search_MXE[get("N") > 1]
    introns_search_MXE_pos <- introns_search_MXE[get("strand") == "+"]
    setorderv(introns_search_MXE_pos,
        c("seqnames", "intron_start", "intron_end"))
    introns_search_MXE_neg <- introns_search_MXE[get("strand") == "-"]
    setorderv(introns_search_MXE_neg,
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_MXE <- rbindlist(
        list(introns_search_MXE_pos, introns_search_MXE_neg))
    introns_search_MXE <- introns_search_MXE[,
        c("skip_coord", "gene_id", "Event", "transcript_id",
            "transcript_name", "intron_number")]
    setnames(introns_search_MXE, old = "Event", new = "Event1")

    introns_search_MXE2 <- introns.skipcoord[,
        c("skip_coord_2", "gene_id", "Event",
        "transcript_id", "transcript_name")]
    setnames(introns_search_MXE2, old = c("skip_coord_2", "Event"),
        new = c("skip_coord", "Event2"))

    introns_search_MXE[introns_search_MXE2,
        on = c("gene_id", "transcript_id", "transcript_name", "skip_coord"),
        c("Event2") := get("i.Event2")]

    introns_search_MXE <- unique(introns_search_MXE,
        by = c("Event1", "Event2"))

    if (nrow(introns_search_MXE) > 0) {
        introns_found_MXE <- introns_search_MXE[,
            {
                edge1 <- rep(seq_len(.N), (.N:1) - 1L)
                i <- 2L:(.N * (.N - 1L) / 2L + 1L)
                o <- cumsum(c(0, (.N - 2L):1))
                edge2 <- i - o[edge1]
                .(
                    gene_id = get("gene_id")[edge1],
                    gene_id_b = get("gene_id")[edge2],
                    Event1a = get("Event1")[edge1],
                    Event1b = get("Event1")[edge2],
                    Event2a = get("Event2")[edge1],
                    Event2b = get("Event2")[edge2],
                    transcript_id_a = get("transcript_id")[edge1],
                    transcript_id_b = get("transcript_id")[edge2],
                    transcript_name_a = get("transcript_name")[edge1],
                    transcript_name_b = get("transcript_name")[edge2],
                    intron_number_a = get("intron_number")[edge1],
                    intron_number_b = get("intron_number")[edge2]
                )
            }, by = "skip_coord"
        ]
        # Make sure to exclude A3SS / A5SS events:
        introns_found_MXE <- introns_found_MXE[get("Event1a") != get("Event1b")]
        introns_found_MXE <- introns_found_MXE[get("Event2a") != get("Event2b")]

        setorderv(introns_found_MXE, c("gene_id", "transcript_name_a"))
        introns_found_MXE[, c("EventName") := paste0(
                "MXE:", get("transcript_name_a"), "-exon",
                (1 + get("intron_number_a")), ";",
                get("transcript_name_b"), "-exon",
                (1 + get("intron_number_b")))
        ]
        introns_found_MXE[, c("EventID") := paste0("MXE#", seq_len(.N))]
        setnames(introns_found_MXE,
            old = "skip_coord", new = "EventRegion")
        introns_found_MXE[, c("EventType") := "MXE"]
        introns_found_MXE <- introns_found_MXE[,
            c("EventType", "EventID", "EventName", "Event1a", "Event1b",
                "Event2a", "Event2b", "gene_id", "gene_id_b", "EventRegion",
                "transcript_id_a", "transcript_name_a", "intron_number_a",
                "transcript_id_b", "transcript_name_b", "intron_number_b")]
        introns_found_MXE <- unique(introns_found_MXE,
            by = c("Event1a", "Event1b", "Event2a", "Event2b"))
    } else {
        introns_found_MXE <- c()
    }
    return(introns_found_MXE)
}

# Generate a list of SE
.gen_splice_SE <- function(introns.skipcoord, candidate.introns) {
    introns.skippedJn <- introns.skipcoord[
        get("skip_coord") %in% get("Event"),
        c("gene_id", "gene_name", "skip_coord")]
    introns.skippedJn <- unique(introns.skippedJn)

    # Get list of skip coords
    introns_skip_SE <- introns.skippedJn[, "skip_coord"]
    introns_search_SE <- candidate.introns[,
        c("gene_id", "Event", "transcript_id",
            "transcript_name", "intron_number")]
    setnames(introns_search_SE,
        old = c("Event", "transcript_id", "transcript_name", "intron_number"),
        new = c("skip_coord", "skip_transcript_id", "skip_transcript_name",
            "skip_intron_number"))
    introns_skip_SE[introns_search_SE, on = "skip_coord",
        c("gene_id_b", "skip_transcript_id",
            "skip_transcript_name", "skip_intron_number") :=
            list(get("i.gene_id"), get("i.skip_transcript_id"),
                get("i.skip_transcript_name"), get("i.skip_intron_number"))
    ]
    introns_skip_SE <- unique(introns_skip_SE,
        by = c("gene_id_b", "skip_coord"))

    # Get junction EiEi+1 and annotate "skip" junction as EiEi+2
    introns_found_SE <- introns.skipcoord[,
        c("skip_coord", "gene_id", "Event", "transcript_id",
            "transcript_name", "intron_number")]
    setnames(introns_found_SE,
        old = c("transcript_id", "intron_number", "transcript_name"),
        new = c("inc_transcript_id", "inc_intron_number",
            "inc_transcript_name"))
    introns_found_SE[introns_skip_SE, on = "skip_coord",
        c("skip_transcript_id", "skip_transcript_name",
            "skip_intron_number", "gene_id_b") :=
            list(get("i.skip_transcript_id"), get("i.skip_transcript_name"),
                get("i.skip_intron_number"), get("i.gene_id_b"))
    ]
    introns_found_SE <- na.omit(introns_found_SE)

    # Same, but "skip" is now Ei-1Ei+1
    introns_found_SE2 <- introns.skipcoord[,
        c("skip_coord_2", "gene_id", "Event",
            "transcript_id", "transcript_name")]
    setnames(introns_found_SE2, old = c("skip_coord_2", "transcript_id"),
        new = c("skip_coord", "inc_transcript_id"))
    introns_found_SE2[introns_skip_SE, on = "skip_coord",
        c("skip_transcript_id",
            "skip_transcript_name", "skip_intron_number") :=
            list(get("i.skip_transcript_id"),
                get("i.skip_transcript_name"), get("i.skip_intron_number"))
    ]
    introns_found_SE2 <- na.omit(introns_found_SE2)

    # Marry the two
    introns_found_SE[introns_found_SE2,
        on = c("skip_coord", "inc_transcript_id"),
        c("inc_coord_downst") := get("i.Event")]
    introns_found_SE <- na.omit(introns_found_SE)

    setnames(introns_found_SE,
        old = c("Event", "inc_coord_downst", "skip_coord"),
        new = c("Event1a", "Event2a", "Event1b")
    )
    introns_found_SE <- unique(introns_found_SE,
        by = c("Event1a", "Event2a", "Event1b"))

    setorderv(introns_found_SE, c("gene_id", "inc_transcript_name"))
    introns_found_SE[, c("EventName") := paste0(
        "SE:", get("inc_transcript_name"), "-exon",
        (1 + get("inc_intron_number")), ";",
        get("skip_transcript_name"), "-int", get("skip_intron_number"))]
    introns_found_SE[, c("EventID") := paste0("SE#", seq_len(.N))]
    introns_found_SE[, c("EventType") := "SE"]
    introns_found_SE[, c("Event2b") := NA]
    introns_found_SE[, c("EventRegion") := get("Event1b")]
    introns_found_SE <- introns_found_SE[,
        c("EventType", "EventID", "EventName", "Event1a",
            "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "inc_transcript_id", "inc_transcript_name", "inc_intron_number",
            "skip_transcript_id", "skip_transcript_name", "skip_intron_number")]
    setnames(introns_found_SE,
        old = c("inc_transcript_id", "inc_transcript_name",
            "skip_transcript_id", "skip_transcript_name"),
        new = c("transcript_id_a", "transcript_name_a",
            "transcript_id_b", "transcript_name_b"))
    setnames(introns_found_SE,
        new = c("intron_number_a", "intron_number_b"),
        old = c("inc_intron_number", "skip_intron_number"))
    return(introns_found_SE)
}

# Generate a list of AFE
.gen_splice_AFE <- function(candidate.introns, introns_found_A5SS) {
    introns_search_AFE <- candidate.introns[get("intron_number") == 1]
    introns_search_AFE_pos <- introns_search_AFE[get("strand") == "+"]
    setorderv(introns_search_AFE_pos,
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_AFE_pos <- introns_search_AFE_pos[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_AFE_pos, old = "intron_end", new = "intron_coord")

    introns_search_AFE_neg <- introns_search_AFE[get("strand") == "-"]
    setorderv(introns_search_AFE_neg,
        c("seqnames", "intron_start", "intron_end"))
    introns_search_AFE_neg <- introns_search_AFE_neg[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_AFE_neg, old = "intron_start", new = "intron_coord")

    introns_search_AFE <- rbindlist(
        list(introns_search_AFE_pos, introns_search_AFE_neg))
    introns_search_AFE <- unique(introns_search_AFE, by = "Event")

    introns_found_AFE <- introns_search_AFE[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2]
            )
        }, by = c("seqnames", "intron_coord")
    ]
    introns_found_AFE <- introns_found_AFE[!is.na(get("gene_id"))]
    introns_found_AFE <- unique(introns_found_AFE, by = c("Event1a", "Event1b"))

    setorderv(introns_found_AFE, c("gene_id", "transcript_name_a"))
    introns_found_AFE <- introns_found_AFE[,
        c("EventName") := paste0(
            "AFE:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b")
        )]
    introns_found_AFE <- introns_found_AFE[, c("EventType") := "AFE"]
    introns_found_AFE <- introns_found_AFE[, c("EventRegion") := get("Event1b")]

    introns_found_AFE <- introns_found_AFE[!introns_found_A5SS,
        on = c("Event1a", "Event1b")]

    introns_found_AFE[, c("EventID") := paste0("AFE#", seq_len(.N))]
    introns_found_AFE <- introns_found_AFE[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_AFE)
}

# Generate a list of ALE
.gen_splice_ALE <- function(candidate.introns, introns_found_A3SS) {
    introns_search_ALE <- candidate.introns[candidate.introns[,
        .I[get("intron_number") == max(get("intron_number"))],
        by = "transcript_id"]$V1]
    introns_search_ALE_pos <- introns_search_ALE[get("strand") == "+"]
    setorderv(introns_search_ALE_pos,
        c("seqnames", "intron_start", "intron_end"))
    introns_search_ALE_pos <- introns_search_ALE_pos[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_ALE_pos,
        old = "intron_start", new = "intron_coord")

    introns_search_ALE_neg <- introns_search_ALE[get("strand") == "-"]
    setorderv(introns_search_ALE_neg,
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_ALE_neg <- introns_search_ALE_neg[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_ALE_neg, old = "intron_end", new = "intron_coord")

    introns_search_ALE <- rbindlist(
        list(introns_search_ALE_pos, introns_search_ALE_neg))
    introns_search_ALE <- unique(introns_search_ALE, by = "Event")

    introns_found_ALE <- introns_search_ALE[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2]
            )
        }, by = c("seqnames", "intron_coord")
    ]
    introns_found_ALE <- introns_found_ALE[!is.na(get("gene_id"))]
    introns_found_ALE <- unique(introns_found_ALE, by = c("Event1a", "Event1b"))

    setorderv(introns_found_ALE, c("gene_id", "transcript_name_a"))
    introns_found_ALE <- introns_found_ALE[,
        c("EventName") := paste0("ALE:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b"))]
    introns_found_ALE <- introns_found_ALE[, c("EventType") := "ALE"]
    introns_found_ALE <- introns_found_ALE[, c("EventRegion") := get("Event1b")]

    introns_found_ALE <- introns_found_ALE[!introns_found_A3SS,
        on = c("Event1a", "Event1b")]

    introns_found_ALE[, c("EventID") := paste0("ALE#", seq_len(.N))]
    introns_found_ALE <- introns_found_ALE[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_ALE)
}

# Generate a list of introns with valid exon island coordinates
.gen_splice_ASS_common <- function(candidate.introns) {
    candidate.introns.ASS <- candidate.introns[
        !is.na(get("exon_group_stranded_upstream")) &
            !is.na(get("exon_group_stranded_downstream"))]
    setnames(candidate.introns.ASS,
        old = c("exon_group_stranded_upstream",
            "exon_group_stranded_downstream"),
        new = c("exon_groups_start", "exon_groups_end"))
    return(candidate.introns.ASS)
}

# Generate a list of A5SS
.gen_splice_A5SS <- function(candidate.introns) {
    introns_search_A5SS <- copy(.gen_splice_ASS_common(candidate.introns))
    introns_search_A5SS_pos <- introns_search_A5SS[get("strand") == "+"]

    setorderv(introns_search_A5SS_pos,
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_A5SS_pos <- introns_search_A5SS_pos[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A5SS_pos,
        old = "intron_end", new = "intron_coord")

    introns_search_A5SS_neg <- introns_search_A5SS[get("strand") == "-"]
    setorderv(introns_search_A5SS_neg,
        c("seqnames", "intron_start", "intron_end"))
    introns_search_A5SS_neg <- introns_search_A5SS_neg[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A5SS_neg,
        old = "intron_start", new = "intron_coord")

    introns_search_A5SS <- rbindlist(
        list(introns_search_A5SS_pos, introns_search_A5SS_neg))
    introns_search_A5SS <- unique(introns_search_A5SS, by = "Event")
    introns_found_A5SS <- introns_search_A5SS[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2],
                exon_groups_start_a = get("exon_groups_start")[edge1],
                exon_groups_start_b = get("exon_groups_start")[edge2],
                exon_groups_end_a = get("exon_groups_end")[edge1],
                exon_groups_end_b = get("exon_groups_end")[edge2]
            )
        }, by = "intron_coord"
    ]
    introns_found_A5SS <- introns_found_A5SS[!is.na(get("gene_id"))]
    introns_found_A5SS <- unique(introns_found_A5SS,
        by = c("Event1a", "Event1b"))
    # filter by same exon group starts and ends:
    introns_found_A5SS <- introns_found_A5SS[
        get("exon_groups_start_a") == get("exon_groups_start_b")]
    introns_found_A5SS <- introns_found_A5SS[
        get("exon_groups_end_a") == get("exon_groups_end_b")]
    setorderv(introns_found_A5SS, c("gene_id", "transcript_name_a"))
    introns_found_A5SS <- introns_found_A5SS[,
        c("EventName") := paste0(
            "A5SS:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b"))]
    introns_found_A5SS <- introns_found_A5SS[, c("EventType") := "A5SS"]
    introns_found_A5SS <- introns_found_A5SS[,
        c("EventRegion") := get("Event1b")]

    introns_found_A5SS[, c("EventID") := paste0("A5SS#", seq_len(.N))]
    introns_found_A5SS <- introns_found_A5SS[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_A5SS)
}

# Generate a list of A3SS
.gen_splice_A3SS <- function(candidate.introns) {
    introns_search_A3SS <- copy(
        .gen_splice_ASS_common(candidate.introns))
    introns_search_A3SS_pos <- introns_search_A3SS[get("strand") == "+"]
    setorderv(introns_search_A3SS_pos,
        c("seqnames", "intron_start", "intron_end"))
    introns_search_A3SS_pos <- introns_search_A3SS_pos[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A3SS_pos,
        old = "intron_start", new = "intron_coord")

    introns_search_A3SS_neg <- introns_search_A3SS[get("strand") == "-"]
    setorderv(introns_search_A3SS_neg,
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_A3SS_neg <- introns_search_A3SS_neg[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A3SS_neg,
        old = "intron_end", new = "intron_coord")

    introns_search_A3SS <- rbindlist(
        list(introns_search_A3SS_pos, introns_search_A3SS_neg))
    introns_search_A3SS <- unique(introns_search_A3SS, by = "Event")

    introns_found_A3SS <- introns_search_A3SS[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2],
                exon_groups_start_a = get("exon_groups_start")[edge1],
                exon_groups_start_b = get("exon_groups_start")[edge2],
                exon_groups_end_a = get("exon_groups_end")[edge1],
                exon_groups_end_b = get("exon_groups_end")[edge2]
            )
        }, by = "intron_coord"
    ]
    introns_found_A3SS <- introns_found_A3SS[!is.na(get("gene_id"))]
    introns_found_A3SS <- unique(introns_found_A3SS,
        by = c("Event1a", "Event1b"))
    # filter by same exon group starts and ends:
    introns_found_A3SS <- introns_found_A3SS[
        get("exon_groups_start_a") == get("exon_groups_start_b")]
    introns_found_A3SS <- introns_found_A3SS[
        get("exon_groups_end_a") == get("exon_groups_end_b")]

    setorderv(introns_found_A3SS, c("gene_id", "transcript_name_a"))
    introns_found_A3SS <- introns_found_A3SS[,
        c("EventName") := paste0(
            "A3SS:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b"))
    ]
    introns_found_A3SS <- introns_found_A3SS[, c("EventType") := "A3SS"]
    introns_found_A3SS <- introns_found_A3SS[,
        c("EventRegion") := get("Event1b")]

    introns_found_A3SS[, c("EventID") := paste0("A3SS#", seq_len(.N))]
    introns_found_A3SS <- introns_found_A3SS[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_A3SS)
}

# Generate a list of RI
.gen_splice_RI <- function(candidate.introns, reference_path) {
    Exons <- .grDT(
        read.fst(file.path(reference_path, "fst", "Exons.fst")),
        keep.extra.columns = TRUE
    )
    candidate.RI <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Introns.Dir.fst")))
    candidate.RI[, c("start", "end") :=
        list(get("intron_start"), get("intron_end"))]
    OL <- findOverlaps(
        .grDT(candidate.RI), Exons, type = "within"
    )
    candidate.RI <- candidate.RI[unique(from(OL))]
    
    found_RI <- data.table(EventType = "RI",
        EventID = paste0("RI#", seq_len(nrow(candidate.RI))),
        EventName = "", # assign later
        Event1a = with(candidate.RI,
            paste0(seqnames, ":", start, "-", end, "/", strand)),
        Event1b = with(candidate.RI,
            paste0(seqnames, ":", start, "-", end, "/", strand)),
        Event2a = NA,
        Event2b = NA,
        gene_id = Exons$gene_id[match(candidate.RI$transcript_id,
            Exons$transcript_id)]
    )
    found_RI$gene_id_b <- found_RI$gene_id
    found_RI$EventRegion <- found_RI$Event1a
    found_RI$transcript_id_a <- ""
    found_RI$transcript_name_a <- ""
    found_RI$intron_number_a <- 0
    found_RI$transcript_id_b <- candidate.RI$transcript_id
    found_RI$transcript_name_b <- Exons$gene_id[match(
        candidate.RI$transcript_id, Exons$transcript_id)]
    found_RI$intron_number_b <- as.numeric(tstrsplit(candidate.RI$intron_id,
        split = "Intron")[[2]])

    return(found_RI)
}

# Renames the ASEs based on the ranking order of transcripts
.gen_splice_save <- function(AS_Table, candidate.introns, reference_path) {
    candidate.introns.order <- copy(candidate.introns)
    if (!("transcript_support_level" %in% colnames(candidate.introns))) {
        candidate.introns.order[, c("transcript_support_level") := "NA"]
    }
    candidate.introns.order[,
        c("is_protein_coding") := !is.na(get("protein_id"))]
    candidate.introns.order[, by = "transcript_id",
        c("is_last_intron") :=
            (get("intron_number") == max(get("intron_number")))]

    AS_Table <- .gen_splice_prep_events(AS_Table, candidate.introns.order,
        reference_path)
    AS_Table <- .gen_splice_name_events(AS_Table, reference_path)
}

.gen_splice_prep_events_RI <- function(AS_Table, Exons) {
    if (!("transcript_support_level" %in% names(mcols(Exons)))) {
        Exons$transcript_support_level <- NA
    }
    Exons$transcript_support_level <- tstrsplit(Exons$transcript_support_level,
        split = " ", fixed = TRUE)[[1]]
    RI.ranges <- AS_Table[get("EventType") == "RI"]
    RI.gr <- coord2GR(RI.ranges$Event1b)
    OL <- findOverlaps(RI.gr, Exons, type = "within")
    RI.DT <- data.table(
        EventType = "RI",
        EventID = RI.ranges$EventID[from(OL)],
        Event1a = RI.ranges$Event1a[from(OL)],
        Event2a = NA,
        transcript_id = Exons$transcript_id[to(OL)]
    )
    RI.DT$transcript_support_level <-
        Exons$transcript_support_level[to(OL)]
    RI.DT$is_protein_coding <-
        (Exons$transcript_biotype[to(OL)] == "protein_coding")
    RI.DT$is_last_intron <- FALSE
    RI.DT$in_1a <- Exons$exon_number[to(OL)]
    RI.DT$in_2a <- Exons$exon_number[to(OL)]
    return(RI.DT)
}

.gen_splice_prep_events <- function(AS_Table, candidate.introns.order,
        reference_path) {
    Exons <- .grDT(
        read.fst(file.path(reference_path, "fst", "Exons.fst")),
        keep.extra.columns = TRUE
    )

    AS_Table_search.a <- AS_Table[get("EventType") != "RI",
        c("EventType", "EventID", "Event1a", "Event2a")]
    AS_Table_search.a[, c("Event") := get("Event1a")]
    AS_Table_search.a <- candidate.introns.order[AS_Table_search.a,
        on = "Event",
        c("EventType", "EventID", "Event1a", "Event2a", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "intron_number")]
    setnames(AS_Table_search.a, "intron_number", "in_1a")
    AS_Table_search.a <-
        AS_Table_search.a[get("EventType") != "AFE" | get("in_1a") == 1]
    AS_Table_search.a <-
        AS_Table_search.a[get("EventType") != "ALE" | get("is_last_intron")]
    AS_Table_search.a[, c("Event") := get("Event2a")]
    AS_Table_search.a[is.na(get("Event")), c("Event") := get("Event1a")]
    AS_Table_search.a <- candidate.introns.order[AS_Table_search.a,
        on = c("Event", "transcript_id", "transcript_support_level"),
        c("EventType", "EventID", "Event1a", "Event2a", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "in_1a", "intron_number")]
    AS_Table_search.a <- AS_Table_search.a[!is.na(get("intron_number"))]
    setnames(AS_Table_search.a, "intron_number", "in_2a")

    # Separate search_a for RI
    RI.DT <- .gen_splice_prep_events_RI(AS_Table, Exons)
    AS_Table_search.a <- rbind(AS_Table_search.a, RI.DT)

    AS_Table_search.b <- AS_Table[,
        c("EventType", "EventID", "Event1b", "Event2b")]
    AS_Table_search.b[, c("Event") := get("Event1b")]
    AS_Table_search.b <- candidate.introns.order[AS_Table_search.b,
        on = "Event",
        c("EventType", "EventID", "Event1b", "Event2b", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "intron_number")]
    setnames(AS_Table_search.b, "intron_number", "in_1b")
    AS_Table_search.b <-
        AS_Table_search.b[get("EventType") != "AFE" | get("in_1b") == 1]
    AS_Table_search.b <-
        AS_Table_search.b[get("EventType") != "ALE" | get("is_last_intron")]
    AS_Table_search.b[, c("Event") := get("Event2b")]
    AS_Table_search.b[is.na(get("Event")), c("Event") := get("Event1b")]
    AS_Table_search.b <- candidate.introns.order[AS_Table_search.b,
        on = c("Event", "transcript_id", "transcript_support_level"),
        c("EventType", "EventID", "Event1b", "Event2b", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "in_1b", "intron_number")]
    AS_Table_search.b <- AS_Table_search.b[!is.na(get("intron_number"))]
    setnames(AS_Table_search.b, "intron_number", "in_2b")

    AS_Table_search.a$transcript_name <- Exons$transcript_name[match(
        AS_Table_search.a$transcript_id, Exons$transcript_id)]
    AS_Table_search.b$transcript_name <- Exons$transcript_name[match(
        AS_Table_search.b$transcript_id, Exons$transcript_id)]
    setorderv(AS_Table_search.a,
        c("transcript_support_level", "is_protein_coding", "transcript_name"),
        order = c(1, -1, 1))
    setorderv(AS_Table_search.b,
        c("transcript_support_level", "is_protein_coding", "transcript_name"),
        order = c(1, -1, 1))
    AS_Table.find.a <- unique(AS_Table_search.a, by = "EventID")
    AS_Table.find.a <- AS_Table.find.a[AS_Table[, "EventID"],
        on = "EventID"]
    AS_Table.find.b <- unique(AS_Table_search.b, by = "EventID")
    AS_Table.find.b <- AS_Table.find.b[AS_Table[, "EventID"],
        on = "EventID"]

    AS_Table$transcript_id_a <- AS_Table.find.a$transcript_id
    AS_Table$transcript_name_a <- AS_Table.find.a$transcript_name
    AS_Table$intron_number_a <- as.numeric(AS_Table.find.a$in_1a)
    AS_Table$transcript_id_b <- AS_Table.find.b$transcript_id
    AS_Table$transcript_name_b <- AS_Table.find.b$transcript_name
    AS_Table$intron_number_b <- as.numeric(AS_Table.find.b$in_1b)

    setnames(AS_Table_search.a,
        old = c("Event1a", "Event2a", "in_1a", "in_2a"),
        new = c("Event1", "Event2", "in_1", "in_2"))
    AS_Table_search.a[, c("isoform") := "A"]
    setnames(AS_Table_search.b,
        old = c("Event1b", "Event2b", "in_1b", "in_2b"),
        new = c("Event1", "Event2", "in_1", "in_2"))
    AS_Table_search.b[, c("isoform") := "B"]

    write.fst(as.data.frame(rbind(AS_Table_search.a, AS_Table_search.b)),
        file.path(reference_path, "fst", "Splice.options.fst"))

    return(AS_Table)
}

# Generate names of splice events
.gen_splice_name_events <- function(AS_Table, reference_path) {
    AS_Table[get("EventType") == "MXE",
        c("EventName") := paste0(
            "MXE:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a")) + 1), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b")) + 1))]
    AS_Table[
        get("EventType") == "SE",
        c("EventName") := paste0(
            "SE:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a")) + 1), ";",
            get("transcript_name_b"), "-int",
            as.character(as.numeric(get("intron_number_b"))))]
    AS_Table[
        get("EventType") == "AFE",
        c("EventName") := paste0(
            "AFE:", get("transcript_name_a"), "-exon1;",
            get("transcript_name_b"), "-exon1")]
    AS_Table[
        get("EventType") == "ALE",
        c("EventName") := paste0(
            "ALE:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a")) + 1), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b")) + 1))]
    AS_Table[
        get("EventType") == "A5SS",
        c("EventName") := paste0(
            "A5SS:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a"))), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b"))))]
    AS_Table[
        get("EventType") == "A3SS",
        c("EventName") := paste0(
            "A3SS:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a") + 1)), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b") + 1)))]
    AS_Table[
        get("EventType") == "RI",
        c("EventName") := paste0(
            "RI:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a"))), ";",
            get("transcript_name_b"), "-intron",
            as.character(as.numeric(get("intron_number_b"))))]
    write.fst(as.data.frame(AS_Table),
        file.path(reference_path, "fst", "Splice.fst"))
    return(AS_Table)
}

################################################################################
# Sub

# Generate nucleotide and peptide sequences for ASE
.gen_splice_proteins <- function(reference_path, genome) {
    .log("Translating Alternate Splice Peptides...",
        "message", appendLF = FALSE)

    AS_Table <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Splice.fst"))
    )
    Proteins_Splice <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Proteins.fst"))
    )
    AS_Table.Extended <- copy(AS_Table)

    Proteins_Splice$exon_number <- as.numeric(Proteins_Splice$exon_number)
    # make phase easier for me to understand
    Proteins_Splice[, c("phase") := -get("phase") %% 3]

    AS_Table.Extended <- .gen_splice_proteins_upstream(AS_Table,
        AS_Table.Extended, Proteins_Splice, genome, isoform = "A")
    AS_Table.Extended <- .gen_splice_proteins_upstream(AS_Table,
        AS_Table.Extended, Proteins_Splice, genome, isoform = "B")
    AS_Table.Extended <- .gen_splice_proteins_downstream(AS_Table,
        AS_Table.Extended, Proteins_Splice, genome, isoform = "A")
    AS_Table.Extended <- .gen_splice_proteins_downstream(AS_Table,
        AS_Table.Extended, Proteins_Splice, genome, isoform = "B")
    AS_Table.Extended <- .gen_splice_proteins_casette(AS_Table,
        AS_Table.Extended, Proteins_Splice, genome, isoform = "A")
    AS_Table.Extended <- .gen_splice_proteins_casette(AS_Table,
        AS_Table.Extended, Proteins_Splice, genome, isoform = "B")

    AS_Table.Extended <- .gen_splice_proteins_trim(AS_Table.Extended)
    AS_Table.Extended <- .gen_splice_proteins_translate(AS_Table.Extended)

    AS_Table.Extended[, c("AA_full_A") := ""]
    AS_Table.Extended[!is.na(get("AA_upstr_A")),
        c("AA_full_A") := paste0(get("AA_full_A"), get("AA_upstr_A"))]
    AS_Table.Extended[!is.na(get("AA_casette_A")),
        c("AA_full_A") := paste0(get("AA_full_A"), get("AA_casette_A"))]
    AS_Table.Extended[!is.na(get("AA_downstr_A")),
        c("AA_full_A") := paste0(get("AA_full_A"), get("AA_downstr_A"))]
    AS_Table.Extended[, c("AA_full_B") := ""]
    AS_Table.Extended[!is.na(get("AA_upstr_B")),
        c("AA_full_B") := paste0(get("AA_full_B"), get("AA_upstr_B"))]
    AS_Table.Extended[!is.na(get("AA_casette_B")),
        c("AA_full_B") := paste0(get("AA_full_B"), get("AA_casette_B"))]
    AS_Table.Extended[!is.na(get("AA_downstr_B")),
        c("AA_full_B") := paste0(get("AA_full_B"), get("AA_downstr_B"))]
    write.fst(as.data.frame(AS_Table.Extended),
        file.path(reference_path, "fst", "Splice.Extended.fst"))
    message("done")
}

# Generate upstream nucleotides
.gen_splice_proteins_upstream <- function(AS_Table, AS_Table.Extended,
        Proteins_Splice, genome, isoform = c("A", "B")) {
    # Upstream applicable for MXE, SE, ALE, A3SS
    cols <- c("EventType", "EventID",
        paste0(c("transcript_id_", "intron_number_"), tolower(isoform)))
    if (isoform == "A") {
        Upstream <- AS_Table[get("EventType") %in%
            c("MXE", "SE", "ALE", "A3SS"),
            cols, with = FALSE]
    } else {
        Upstream <- AS_Table[get("EventType") %in%
            c("MXE", "SE", "ALE", "A3SS", "RI"),
            cols, with = FALSE]
    }

    Upstream[, c("transcript_id", "exon_number") :=
        list(get(paste0("transcript_id_", tolower(isoform))),
            get(paste0("intron_number_", tolower(isoform))))]

    Upstream <- Proteins_Splice[Upstream,
        on = c("transcript_id", "exon_number"),
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")
    ] # left_join with Exons
    Upstream_gr <- .grDT(na.omit(Upstream), keep.extra.columns = TRUE)

    Upstream_seq <- getSeq(genome, Upstream_gr)
    Upstream[!is.na(get("seqnames")), c("seq") := as.character(Upstream_seq)]

    cols <- c("EventID", "phase", "seq")
    AS_Table.Extended[Upstream[, cols, with = FALSE], on = "EventID",
        c(paste0(c("phase_upstr_", "DNA_upstr_"), toupper(isoform))) :=
        list(get("phase"), get("seq"))]
    return(AS_Table.Extended)
}

# Generate downstream nucleotides
.gen_splice_proteins_downstream <- function(AS_Table, AS_Table.Extended,
        Proteins_Splice, genome, isoform = c("A", "B")) {
    # Add EventType as exon_number is conditional on this
    cols <- c("EventType", "EventID",
        paste0(c("transcript_id_", "intron_number_"), tolower(isoform)))
    if (isoform == "A") {
        Downstream <- AS_Table[get("EventType") %in%
            c("MXE", "SE", "AFE", "A5SS"),
            cols, with = FALSE]
    } else {
        Downstream <- AS_Table[get("EventType") %in%
            c("MXE", "SE", "AFE", "A5SS", "RI"),
            cols, with = FALSE]
    }

    Downstream[, c("transcript_id", "exon_number") :=
        list(get(paste0("transcript_id_", tolower(isoform))),
            get(paste0("intron_number_", tolower(isoform))))]
    # Modify downstream exon number
    if (toupper(isoform) == "A") {
        Downstream[get("EventType") %in% c("MXE", "SE"),
            c("exon_number") := get("exon_number") + 2]
        Downstream[get("EventType") %in% c("AFE", "A5SS"),
            c("exon_number") := get("exon_number") + 1]
    } else {
        Downstream[get("EventType") %in% c("MXE"),
            c("exon_number") := get("exon_number") + 2]
        Downstream[get("EventType") %in% c("SE", "AFE", "A5SS", "RI"),
            c("exon_number") := get("exon_number") + 1]
    }

    # left_join with Exons
    Downstream <- Proteins_Splice[Downstream,
        on = c("transcript_id", "exon_number"),
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")
    ] # left_join with Exons
    Downstream_gr <- .grDT(na.omit(Downstream), keep.extra.columns = TRUE)
    Downstream_seq <- getSeq(genome, Downstream_gr)
    Downstream[!is.na(get("seqnames")),
        c("seq") := as.character(Downstream_seq)]

    cols <- c("EventID", "phase", "seq")
    AS_Table.Extended[Downstream[, cols, with = FALSE], on = "EventID",
        c(paste0(c("phase_downstr_", "DNA_downstr_"),
            toupper(isoform))) :=
        list(get("phase"), get("seq"))] # translate
    return(AS_Table.Extended)
}

# Generate casette nucleotides
.gen_splice_proteins_casette <- function(AS_Table, AS_Table.Extended,
        Proteins_Splice, genome, isoform = c("A", "B")) {
    cols <- c("EventType", "EventID",
        paste0(c("transcript_id_", "intron_number_"), tolower(isoform)))
    Casette <- AS_Table[, cols, with = FALSE]
    Casette[, c("transcript_id", "exon_number") :=
        list(get(paste0("transcript_id_", tolower(isoform))),
            get(paste0("intron_number_", tolower(isoform))))]
    if (toupper(isoform) == "B") {
        Casette <- Casette[!(get("EventType") %in% c("SE", "RI"))]
        Casette[get("EventType") %in% c("MXE", "ALE", "A3SS"),
            c("exon_number") := get("exon_number") + 1]
    } else {
        Casette[get("EventType") %in% c("MXE", "SE", "ALE", "A3SS"),
            c("exon_number") := get("exon_number") + 1]
    }
    Casette <- Proteins_Splice[Casette,
        on = c("transcript_id", "exon_number"),
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Casette_gr <- .grDT(na.omit(Casette), keep.extra.columns = TRUE)
    Casette_seq <- getSeq(genome, Casette_gr)
    Casette[!is.na(get("seqnames")),
        c("seq") := as.character(Casette_seq)]

    cols <- c("EventID", "phase", "seq")
    AS_Table.Extended[Casette[, cols, with = FALSE], on = "EventID",
        c(paste0(c("phase_casette_", "DNA_casette_"),
        toupper(isoform))) :=
        list(get("phase"), get("seq"))]
    return(AS_Table.Extended)
}

# Internal trim functions

# Trim 5'-nucleotides based on phase
.trim_phase <- function(DNAstr, phase)
    substr(DNAstr, 1 + (3 - phase) %% 3, nchar(DNAstr))

# Trim 3'-nucleotides based on phase
.trim_3 <- function(DNAstr)
    substr(DNAstr, 1, nchar(DNAstr) - (nchar(DNAstr) %% 3))

# Returns the incomplete codon of nucleotides from the 3'-sequence
.transfer_down <- function(DNAstr)
    substr(DNAstr, nchar(DNAstr) - (nchar(DNAstr) %% 3) + 1, nchar(DNAstr))

# Returns the incomplete codon of nucleotides from the 5'-sequence
.transfer_up <- function(DNAstr, phase)
    substr(DNAstr, 1, 3 - phase)


# Trim nucleotide sequences of upstream / casette / downstream as per phase
.gen_splice_proteins_trim <- function(AS_Table.Extended) {

    AS_Table.Extended <- .gen_splice_proteins_trim_5prime(AS_Table.Extended)

    AS_Table.Extended <- .gen_splice_proteins_transfers(AS_Table.Extended)

    AS_Table.Extended <- .gen_splice_proteins_trim_3prime(AS_Table.Extended)

    return(AS_Table.Extended)
}

# Trims 5'-sequences from upstr, casette, or downstream if there is no
# sequences upstream to it
.gen_splice_proteins_trim_5prime <- function(AS_Table.Extended) {
    # Trim 5' upstream
    AS_Table.Extended[!is.na(get("DNA_upstr_A")) &
        get("phase_upstr_A") != 0,
        c("DNA_upstr_A") := .trim_phase(
            get("DNA_upstr_A"), get("phase_upstr_A"))
    ]
    AS_Table.Extended[!is.na(get("DNA_upstr_B")) &
        get("phase_upstr_B") != 0,
        c("DNA_upstr_B") := .trim_phase(
            get("DNA_upstr_B"), get("phase_upstr_B"))
    ]
    # Trim 5' casette if upstream does not exist
    AS_Table.Extended[
        is.na(get("DNA_upstr_A")) & !is.na(get("DNA_casette_A")) &
        get("phase_casette_A") != 0,
        c("DNA_casette_A") := .trim_phase(
            get("DNA_casette_A"), get("phase_casette_A"))
    ]
    AS_Table.Extended[
        is.na(get("DNA_upstr_B")) & !is.na(get("DNA_casette_B")) &
        get("phase_casette_B") != 0,
        c("DNA_casette_B") := .trim_phase(
            get("DNA_casette_B"), get("phase_casette_B"))
    ]
    # Trim 5' downstream if upstream and casette doesn't exist
    AS_Table.Extended[
        is.na(get("DNA_upstr_A")) & is.na(get("DNA_casette_A")) &
        !is.na(get("DNA_downstr_A")) & get("phase_downstr_A") != 0,
        c("DNA_downstr_A") := .trim_phase(
            get("DNA_downstr_A"), get("phase_downstr_A"))
    ]
    AS_Table.Extended[
        is.na(get("DNA_upstr_B")) & is.na(get("DNA_casette_B")) &
        !is.na(get("DNA_downstr_B")) & get("phase_downstr_B") != 0,
        c("DNA_downstr_B") := .trim_phase(
            get("DNA_downstr_B"), get("phase_downstr_B"))
    ]
    return(AS_Table.Extended)
}

# Transfers nucleotides to neighboring element so that all elements have
# in-frame sequences
.gen_splice_proteins_transfers <- function(AS_Table.Extended) {
    # Transfer from upstream to casette if both exist
    AS_Table.Extended[
        !is.na(get("DNA_upstr_A")) & !is.na(get("DNA_casette_A")) &
        get("phase_casette_A") != 0,
        c("DNA_casette_A") := paste0(
            .transfer_down(get("DNA_upstr_A")), get("DNA_casette_A"))
    ]
    AS_Table.Extended[
        !is.na(get("DNA_upstr_B")) & !is.na(get("DNA_casette_B")) &
        get("phase_casette_B") != 0,
        c("DNA_casette_B") := paste0(
            .transfer_down(get("DNA_upstr_B")), get("DNA_casette_B"))
    ]
    # Transfer from upstream to downstream if both exist and casette doesn't
    AS_Table.Extended[
        !is.na(get("DNA_upstr_A")) & !is.na(get("DNA_downstr_A")) &
        is.na(get("DNA_casette_A")) & get("phase_downstr_A") != 0,
        c("DNA_downstr_A") := paste0(
            .transfer_down(get("DNA_upstr_A")), get("DNA_downstr_A"))
    ]
    AS_Table.Extended[
        !is.na(get("DNA_upstr_B")) & !is.na(get("DNA_downstr_B")) &
        is.na(get("DNA_casette_B")) & get("phase_downstr_B") != 0,
        c("DNA_downstr_B") := paste0(
            .transfer_down(get("DNA_upstr_B")), get("DNA_downstr_B"))
    ]
    # Transfer from downstream to casette if both exist
    AS_Table.Extended[
        !is.na(get("DNA_casette_A")) & !is.na(get("DNA_downstr_A")) &
        get("phase_downstr_A") != 0,
        c("DNA_casette_A") := paste0(get("DNA_casette_A"),
            .transfer_up(get("DNA_downstr_A"), get("phase_downstr_A")))
    ]
    AS_Table.Extended[
        !is.na(get("DNA_casette_B")) & !is.na(get("DNA_downstr_B")) &
        get("phase_downstr_B") != 0,
        c("DNA_casette_B") := paste0(get("DNA_casette_B"),
            .transfer_up(get("DNA_downstr_B"), get("phase_downstr_B")))
    ]
    return(AS_Table.Extended)
}

# Trim 3' ends so translation doesnt return error
.gen_splice_proteins_trim_3prime <- function(AS_Table.Extended) {
    AS_Table.Extended[!is.na(get("DNA_upstr_A")),
        c("DNA_upstr_A") := .trim_3(get("DNA_upstr_A"))]
    AS_Table.Extended[!is.na(get("DNA_casette_A")),
        c("DNA_casette_A") := .trim_3(get("DNA_casette_A"))]
    AS_Table.Extended[!is.na(get("DNA_downstr_A")),
        c("DNA_downstr_A") := .trim_3(get("DNA_downstr_A"))]
    AS_Table.Extended[!is.na(get("DNA_upstr_B")),
        c("DNA_upstr_B") := .trim_3(get("DNA_upstr_B"))]
    AS_Table.Extended[!is.na(get("DNA_casette_B")),
        c("DNA_casette_B") := .trim_3(get("DNA_casette_B"))]
    AS_Table.Extended[!is.na(get("DNA_downstr_B")),
        c("DNA_downstr_B") := .trim_3(get("DNA_downstr_B"))]
    return(AS_Table.Extended)
}

# Translate upstream, casette and downstream sequences
.gen_splice_proteins_translate <- function(AS_Table.Extended) {
    DNAseq <- AS_Table.Extended[nchar(get("DNA_upstr_A")) > 0]$DNA_upstr_A
    AAseq <- as.character(Biostrings::translate(as(DNAseq, "DNAStringSet")))
    AS_Table.Extended[nchar(get("DNA_upstr_A")) > 0,
        c("AA_upstr_A") := AAseq]

    DNAseq <- AS_Table.Extended[nchar(get("DNA_casette_A")) > 0]$DNA_casette_A
    AAseq <- as.character(Biostrings::translate(as(DNAseq, "DNAStringSet")))
    AS_Table.Extended[nchar(get("DNA_casette_A")) > 0,
        c("AA_casette_A") := AAseq]

    DNAseq <- AS_Table.Extended[nchar(get("DNA_downstr_A")) > 0]$DNA_downstr_A
    AAseq <- as.character(Biostrings::translate(as(DNAseq, "DNAStringSet")))
    AS_Table.Extended[nchar(get("DNA_downstr_A")) > 0,
        c("AA_downstr_A") := AAseq]

    DNAseq <- AS_Table.Extended[nchar(get("DNA_upstr_B")) > 0]$DNA_upstr_B
    AAseq <- as.character(Biostrings::translate(as(DNAseq, "DNAStringSet")))
    AS_Table.Extended[nchar(get("DNA_upstr_B")) > 0,
        c("AA_upstr_B") := AAseq]

    DNAseq <- AS_Table.Extended[nchar(get("DNA_casette_B")) > 0]$DNA_casette_B
    AAseq <- as.character(Biostrings::translate(as(DNAseq, "DNAStringSet")))
    AS_Table.Extended[nchar(get("DNA_casette_B")) > 0,
        c("AA_casette_B") := AAseq]

    DNAseq <- AS_Table.Extended[nchar(get("DNA_downstr_B")) > 0]$DNA_downstr_B
    AAseq <- as.character(Biostrings::translate(as(DNAseq, "DNAStringSet")))
    AS_Table.Extended[nchar(get("DNA_downstr_B")) > 0,
        c("AA_downstr_B") := AAseq]

    return(AS_Table.Extended)
}
