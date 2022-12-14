#' STAR wrappers for building reference for STAR, and aligning RNA-sequencing
#'
#' These functions run the STAR aligner to build a STAR genome reference,
#' calculate mappability exclusion regions using STAR, and align one or more
#' FASTQ files (single or paired) to the generated genome. These functions only
#' work on Linux-based systems with STAR installed. STAR must be
#' accessible via `$PATH`. See details and examples
#'
#' @details
#' **Pre-requisites**
#'
#' `STAR_buildRef` requires [getResources] to be run to fetch the
#' required genome and gene annotation files.
#'
#' `STAR_mappability`, `STAR_alignExperiment` and `STAR_alignReads` requires a
#' `STAR` genome, which can be built using `STAR_buildRef`
#'
#' **Function Description**
#'
#' For `STAR_buildRef`: this function
#'   will create a `STAR` genome reference in the `STAR` subdirectory in the
#'   path given by `reference_path`. Optionally, it will run `STAR_mappability`
#'   if `also_generate_mappability` is set to `TRUE`
#'
#' For `STAR_mappability`: this function will first
#'   will run [generateSyntheticReads], then use the given `STAR` genome to 
#'   align the synthetic reads using `STAR`. The aligned BAM file will then be
#'   processed using [calculateMappability] to calculate the
#'   lowly-mappable genomic regions,
#'   producing the `MappabilityExclusion.bed.gz` output file.
#'
#' For `STAR_alignReads`: aligns a single or pair of FASTQ files to the given
#'   `STAR` genome using the `STAR` aligner.
#'
#' For `STAR_alignExperiment`: aligns a set of FASTQ or paired FASTQ files
#'   using the given
#'   `STAR` genome using the `STAR` aligner.
#'   A data.frame specifying sample names and corresponding FASTQ files are
#'   required
#'
#' @param reference_path The path to the reference.
#'    [getResources] must first be run using this path
#'    as its `reference_path`
#' @param STAR_ref_path (Default - the "STAR" subdirectory under
#'    \code{reference_path}) The directory containing the STAR reference to be
#'    used or to contain the newly-generated STAR reference
#' @param also_generate_mappability Whether \code{STAR_buildRef()} also
#'   calculates Mappability Exclusion regions.
#' @param map_depth_threshold (Default 4) The depth of mapped reads
#'   threshold at or below which Mappability exclusion regions are defined. See
#'   [Mappability-methods].
#'   Ignored if \code{also_generate_mappability = FALSE}
#' @param sjdbOverhang (Default = 149) A STAR setting indicating the length of
#'   the donor / acceptor sequence on each side of the junctions. Ideally equal
#'   to (mate_length - 1). As the most common read length is 150, the default
#'   of this function is 149. See the STAR aligner manual for details.
#' @param n_threads The number of threads to run the STAR aligner.
#' @param additional_args A character vector of additional arguments to be
#'   parsed into STAR. See examples below.
#' @param Experiment A two or three-column data frame with the columns denoting
#'   sample names, forward-FASTQ and reverse-FASTQ files. This can be
#'   conveniently generated using [findFASTQ]
#' @param BAM_output_path The path under which STAR outputs the aligned BAM
#'   files. In `STAR_alignExperiment()`, STAR will output aligned
#'   BAMS inside subdirectories of this folder, named by sample names. In
#'   `STAR_alignReads()`, STAR will output directly into this path.
#' @param trim_adaptor The sequence of the Illumina adaptor to trim via STAR's
#'   \code{--clip3pAdapterSeq} option
#' @param two_pass Whether to use two-pass mapping. In
#'   \code{STAR_alignExperiment()}, STAR will first align every sample
#'   and generate a list of splice junctions but not BAM files. The junctions
#'   are then given to STAR to generate a temporary genome (contained within
#'   \code{_STARgenome}) subdirectory within that of the first sample), using
#'   these junctions to improve novel junction detection. In
#'   \code{STAR_alignReads()}, STAR will run \code{--twopassMode Basic}
#' @param fastq_1,fastq_2 In STAR_alignReads: character vectors giving the
#'   path(s) of one or more FASTQ (or FASTA) files to be aligned.
#'   If single reads are to be aligned, omit \code{fastq_2}
#' @param memory_mode The parameter to be parsed to \code{--genomeLoad}; either
#'   \code{NoSharedMemory} or \code{LoadAndKeep} are used.
#' @param overwrite (default `FALSE`) If BAM file(s) already exist from a
#'   previous run, whether these would be overwritten.
#' @param ... Additional arguments to be parsed into
#'   \code{generateSyntheticReads()}. See \link{Mappability-methods}.
#' @return None. STAR will output files into the given output directories.
#' @examples
#' # 0) Check that STAR is installed and compatible with SpliceWiz
#'
#' STAR_version()
#' \dontrun{
#'
#' # The below workflow illustrates
#' # 1) Getting the reference resource
#' # 2) Building the STAR Reference, including Mappability Exclusion calculation
#' # 3) Building the SpliceWiz Reference, using the Mappability Exclusion file
#' # 4) Aligning (a) one or (b) multiple raw sequencing samples.
#'
#'
#' # 1) Reference generation from Ensembl's FTP links
#'
#' FTP <- "ftp://ftp.ensembl.org/pub/release-94/"
#'
#' getResources(
#'     reference_path = "Reference_FTP",
#'     fasta = paste0(FTP, "fasta/homo_sapiens/dna/",
#'         "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
#'     gtf = paste0(FTP, "gtf/homo_sapiens/",
#'         "Homo_sapiens.GRCh38.94.chr.gtf.gz")
#' )
#'
#' # 2) Generates STAR genome within the SpliceWiz reference. Also generates
#' # mappability exclusion gzipped BED file inside the "Mappability/" sub-folder
#'
#' STAR_buildRef(
#'     reference_path = "Reference_FTP",
#'     n_threads = 8,
#'     also_generate_mappability = TRUE
#' )
#'
#' # 2 alt) Generates STAR genome of the example SpliceWiz genome.
#' #     This demonstrates using custom STAR parameters, as the example 
#' #     SpliceWiz genome is ~100k in length, 
#' #     so --genomeSAindexNbases needs to be
#' #     adjusted to be min(14, log2(GenomeLength)/2 - 1)
#'
#' getResources(
#'     reference_path = "Reference_chrZ",
#'     fasta = chrZ_genome(),
#'     gtf = chrZ_gtf()
#' )
#'
#' STAR_buildRef(
#'     reference_path = "Reference_chrZ",
#'     n_threads = 8,
#'     additional_args = c("--genomeSAindexNbases", "7"),
#'     also_generate_mappability = TRUE
#' )
#'
#' # 3) Build SpliceWiz reference using the newly-generated 
#' #    Mappability exclusions
#'
#' #' NB: also specifies to use the hg38 nonPolyA resource
#'
#' buildRef(reference_path = "Reference_FTP", genome_type = "hg38")
#'
#' # 4a) Align a single sample using the STAR reference
#'
#' STAR_alignReads(
#'     STAR_ref_path = file.path("Reference_FTP", "STAR"),
#'     BAM_output_path = "./bams/sample1",
#'     fastq_1 = "sample1_1.fastq", fastq_2 = "sample1_2.fastq",
#'     n_threads = 8
#' )
#'
#' # 4b) Align multiple samples, using two-pass alignment
#'
#' Experiment <- data.frame(
#'     sample = c("sample_A", "sample_B"),
#'     forward = file.path("raw_data", c("sample_A", "sample_B"),
#'         c("sample_A_1.fastq", "sample_B_1.fastq")),
#'     reverse = file.path("raw_data", c("sample_A", "sample_B"),
#'         c("sample_A_2.fastq", "sample_B_2.fastq"))
#' )
#'
#' STAR_alignExperiment(
#'     Experiment = Experiment,
#'     STAR_ref_path = file.path("Reference_FTP", "STAR"),
#'     BAM_output_path = "./bams",
#'     two_pass = TRUE,
#'     n_threads = 8
#' )
#' }
#' @name STAR-methods
#' @aliases
#' STAR_buildRef STAR_alignExperiment STAR_alignReads
#' @seealso
#' [Build-Reference-methods] [findSamples] [Mappability-methods]\cr\cr
#' [The latest STAR documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
#' @md
NULL

#' @describeIn STAR-methods Checks whether STAR is installed, and its version
#' @export
STAR_version <- function() .validate_STAR_version(type = "message")

#' @describeIn STAR-methods Creates a STAR genome reference.
#' @export
STAR_buildRef <- function(reference_path,
        STAR_ref_path = file.path(reference_path, "STAR"),
        also_generate_mappability = TRUE,
        map_depth_threshold = 4,
        sjdbOverhang = 149,
        n_threads = 4,
        additional_args = NULL,
        ...
) {
    .validate_reference_resource(reference_path)
    .validate_STAR_version()
    .validate_path(STAR_ref_path)
    # Unzip reference files
    genome.fa <- .STAR_get_FASTA(reference_path)
    transcripts.gtf <- .STAR_get_GTF(reference_path)
    # Build STAR using defaults
    .log(paste("Building STAR genome from", reference_path), type = "message")

    args <- NULL
    if (!("--runMode" %in% additional_args)) args <- c(
        "--runMode", "genomeGenerate")

    if (!("--genomeDir" %in% additional_args)) args <- c(args,
        "--genomeDir", STAR_ref_path)

    if (!("--genomeFastaFiles" %in% additional_args)) args <- c(args,
        "--genomeFastaFiles", genome.fa)

    if (!("--sjdbGTFfile" %in% additional_args)) args <- c(args,
        "--sjdbGTFfile", transcripts.gtf)

    if (!("--sjdbOverhang" %in% additional_args)) args <- c(args,
        "--sjdbOverhang", sjdbOverhang)

    if (!("--runThreadN" %in% additional_args)) args <- c(args,
        "--runThreadN", .validate_threads(n_threads, as_BPPARAM = FALSE))

    if (!is.null(additional_args) && all(is.character(additional_args))) {
        args <- c(args, additional_args)
    }
    system2(command = "STAR", args = args)

    if (also_generate_mappability) {
        STAR_mappability(
            reference_path = reference_path,
            STAR_ref_path = STAR_ref_path,
            map_depth_threshold = map_depth_threshold,
            n_threads = n_threads,
            ...
        )
    }

    # Clean up
    .STAR_clean_temp_FASTA_GTF(reference_path)
}

#' @describeIn STAR-methods Calculates lowly-mappable genomic regions using STAR
#' @export
STAR_mappability <- function(
        reference_path,
        STAR_ref_path = file.path(reference_path, "STAR"),
        map_depth_threshold = 4,
        n_threads = 4,
        ...
) {
    .validate_reference_resource(reference_path)
    .validate_STAR_version()
    STAR_ref_path <- .validate_STAR_reference(STAR_ref_path)
    mappability_reads_fasta <- file.path(
        reference_path, "Mappability", "Reads.fa")
    generateSyntheticReads(reference_path, ...)

    .log(paste("Aligning genome fragments back to the genome, from:",
        mappability_reads_fasta), type = "message")
    aligned_bam <- file.path(reference_path, "Mappability",
        "Aligned.out.bam")
    STAR_alignReads(
        fastq_1 = mappability_reads_fasta,
        fastq_2 = NULL,
        STAR_ref_path = STAR_ref_path,
        BAM_output_path = dirname(aligned_bam),
        n_threads = n_threads,
        trim_adaptor = "",
        additional_args = c(
            "--outSAMstrandField", "None",
            "--outSAMattributes", "None"
        )
    )
    if (file.exists(aligned_bam)) {
        # Cleaan up fasta
        if (file.exists(mappability_reads_fasta))
            file.remove(mappability_reads_fasta)

        .log(paste("Calculating Mappability from:", aligned_bam),
            type = "message")
        calculateMappability(
            reference_path = reference_path,
            aligned_bam = aligned_bam,
            threshold = map_depth_threshold,
            n_threads = n_threads
        )
    } else {
        .log("STAR failed to align mappability reads", "warning")
    }
    if (file.exists(file.path(reference_path, "Mappability",
            "MappabilityExclusion.bed.gz"))) {
        message("Mappability Exclusion calculations complete")
        # Clean up BAM
        if (file.exists(aligned_bam)) file.remove(aligned_bam)
    } else {
        .log("Mappability Exclusion calculations not performed", "warning")
    }
}

#' @describeIn STAR-methods Aligns multiple sets of FASTQ files, belonging to
#'   multiple samples
#' @export
STAR_alignExperiment <- function(
    Experiment, STAR_ref_path, BAM_output_path,
    trim_adaptor = "AGATCGGAAG", two_pass = FALSE, n_threads = 4,
    overwrite = FALSE
) {

    .validate_STAR_version()
    STAR_ref_path <- .validate_STAR_reference(STAR_ref_path)
    BAM_output_path <- .validate_path(BAM_output_path)

    # Dissect Experiment:
    if (ncol(Experiment) < 2 || ncol(Experiment) > 3) {
        .log(paste("Experiment must be a 2- or 3- column data frame,",
            "with the columns denoting sample name, fastq file (forward),",
            "and (optionally) fastq file (reverse)"))
    } else if (ncol(Experiment) == 2) {
        colnames(Experiment) <- c("sample", "forward")
        fastq_1 <- Experiment[, "forward"]
        fastq_2 <- NULL
        .validate_STAR_fastq_samples(fastq_1)
        paired <- FALSE
    } else if (ncol(Experiment) == 3) {
        colnames(Experiment) <- c("sample", "forward", "reverse")
        fastq_1 <- Experiment[, "forward"]
        fastq_2 <- Experiment[, "reverse"]
        .validate_STAR_fastq_samples(fastq_1, fastq_2)
        paired <- TRUE
    }
    gzipped <- all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if (is_valid(trim_adaptor)) .validate_STAR_trim_sequence(trim_adaptor)

    samples <- unique(Experiment[, "sample"])
    SJ.files <- NULL
    two_pass_genome <- NULL
    loaded_ref <- NULL
    for (pass in seq_len(ifelse(two_pass, 2, 1))) {
        if (two_pass && pass == 1) message("STAR - first pass")
        if (two_pass && pass == 2) message("STAR - second pass")
        if (pass == 1) {
            ref <- STAR_ref_path
            system2(command = "STAR", args = c(
                "--genomeLoad", "LoadAndExit", "--genomeDir", ref,
                "--outFileNamePrefix", tempdir()
            ))
            loaded_ref <- ref
        }
        for (i in seq_len(length(samples))) {
            if (pass == 2 && !is.null(two_pass_genome) && is.null(loaded_ref)) {
                ref <- two_pass_genome
                system2(command = "STAR", args = c(
                    "--genomeLoad", "LoadAndExit", "--genomeDir", ref,
                    "--outFileNamePrefix", tempdir()
                ))
                loaded_ref <- ref
            }

            sample <- samples[i]
            Expr_sample <- Experiment[Experiment[, "sample"] == sample, ]
            if (!paired) {
                fastq_1 <- Expr_sample[, "forward"]
                fastq_2 <- NULL
            } else {
                fastq_1 <- Expr_sample[, "forward"]
                fastq_2 <- Expr_sample[, "reverse"]
            }
            memory_mode <- "LoadAndKeep"
            if (two_pass && pass == 1) {
                additional_args <- c("--outSAMtype", "None")
            } else if (two_pass && pass == 2 && !is.null(SJ.files)) {
                additional_args <- c("--sjdbFileChrStartEnd",
                    paste(SJ.files$path, collapse = " "),
                    "--sjdbInsertSave", "All"
                )
                two_pass_genome <- file.path(BAM_output_path, sample,
                    "_STARgenome")
                SJ.files <- NULL
                memory_mode <- "NoSharedMemory"
            } else {
                additional_args <- NULL
            }

            .log(paste("Aligning", sample, "using STAR"), "message")
            STAR_alignReads(
                STAR_ref_path = ref,
                BAM_output_path = file.path(BAM_output_path, sample),
                fastq_1 = fastq_1, fastq_2 = fastq_2,
                trim_adaptor = trim_adaptor,
                memory_mode = memory_mode,
                additional_args = additional_args,
                n_threads = n_threads,
                overwrite = overwrite
            )

        } # end of FOR loop

        if (two_pass && pass == 1) {
            SJ.files <- findSamples(BAM_output_path, suffix = ".out.tab")
            if (nrow(SJ.files) == 0) {
                .log(paste("In STAR two-pass,",
                    "no SJ.out.tab files were found"))
            }
        }
        .log(paste("Unloading STAR reference:", loaded_ref), "message")
        system2(command = "STAR", args = c(
            "--genomeLoad", "Remove", "--genomeDir", loaded_ref,
                "--outFileNamePrefix", tempdir()
        ))
        loaded_ref <- NULL
    }
}

#' @describeIn STAR-methods Aligns a single sample (with single or paired FASTQ
#'   or FASTA files)
#' @export
STAR_alignReads <- function(
        fastq_1 = c("./sample_1.fastq"), fastq_2 = NULL,
        STAR_ref_path, BAM_output_path,
        two_pass = FALSE,
        trim_adaptor = "AGATCGGAAG",
        memory_mode = "NoSharedMemory",
        additional_args = NULL,
        n_threads = 4,
        overwrite = FALSE
) {
    expectedBAM <- file.path(BAM_output_path, "Aligned.out.bam")
    if(!overwrite) {
        if(file.exists(expectedBAM)) {
            .log(paste(
                expectedBAM, 
                "already exists. Set overwrite = TRUE to overwrite"
                ), "warning"
            )
            message("") # for \n
            return()
        }
    } else {
        .log(paste(
            expectedBAM, "found, overwriting..."
            ), "warning"
        )
        unlink(expectedBAM)
    }

    .validate_STAR_version()
    STAR_ref_path <- .validate_STAR_reference(STAR_ref_path)
    .validate_STAR_fastq_samples(fastq_1, fastq_2)

    paired <- (length(fastq_1) == length(fastq_2))
    gzipped <- all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if (is_valid(trim_adaptor)) .validate_STAR_trim_sequence(trim_adaptor)

    BAM_output_path <- .validate_path(BAM_output_path)
    # Load STAR reference

    # Remove duplication:
    args <- NULL
    if (!("--genomeLoad" %in% additional_args)) 
        args <- c("--genomeLoad", memory_mode)
    if (!("--runThreadN" %in% additional_args)) 
        args <- c(args,
            "--runThreadN", .validate_threads(n_threads, as_BPPARAM = FALSE))
    if (!("--genomeDir" %in% additional_args)) 
        args <- c(args, "--genomeDir", STAR_ref_path)
    if (!("--outFileNamePrefix" %in% additional_args))
        args <- c(args, "--outFileNamePrefix",
            paste0(BAM_output_path, "/"))

    if (!("--outStd" %in% additional_args)) args <- c(args, "--outStd", "Log")
    if (!("--outBAMcompression" %in% additional_args)) args <- c(args, 
        "--outBAMcompression", "6")

    if (!("--outSAMstrandField" %in% additional_args)) args <- c(args, 
        "--outSAMstrandField", "intronMotif")

    if (!("--outSAMunmapped" %in% additional_args)) 
        args <- c(args, "--outSAMunmapped", "None")

    if (!("--outFilterMultimapNmax" %in% additional_args))
        args <- c(args, "--outFilterMultimapNmax", "1")

    if (!("--outSAMtype" %in% additional_args))
        args <- c(args, "--outSAMtype", "BAM", "Unsorted")

    if (two_pass) args <- c(args, "--twopassMode", "Basic")

    args <- c(args, "--readFilesIn", paste(fastq_1, collapse = ","))

    if (paired) args <- c(args, paste(fastq_2, collapse = ","))
    if (gzipped) args <- c(args, "--readFilesCommand", shQuote("gzip -dc"))
    if (is_valid(trim_adaptor)) 
        args <- c(args, "--clip3pAdapterSeq", trim_adaptor)

    if (!is.null(additional_args) && all(is.character(additional_args)))
        args <- c(args, additional_args)

    system2(command = "STAR", args = args)
}

.validate_STAR_version <- function(type = "error") {
    if (Sys.info()["sysname"] != "Linux") {
        .log("STAR wrappers are only supported on Linux", type = type)
    }
    star_version <- NULL
    tryCatch({
        star_version <- system2("STAR", "--version", stdout = TRUE)
    }, error = function(e) {
        star_version <- NULL
    })
    if (is.null(star_version)) {
        .log("STAR is not installed", type = type)
    }
    if (!is.null(star_version) && star_version < "2.5.0") {
        .log(paste("STAR version < 2.5.0 is not supported;",
            "current version:", star_version), type = type)
    }
    if (!is.null(star_version) && star_version >= "2.5.0") {
        .log(paste("STAR version", star_version), type = "message")
    }
}

.validate_STAR_reference <- function(STAR_ref_path) {
    if (!file.exists(file.path(STAR_ref_path, "genomeParameters.txt")))
        .log(paste(
            STAR_ref_path, "does not appear to be a valid STAR reference"))
    return(normalizePath(STAR_ref_path))
}

# Creates a temporary FASTA file from locally-stored TwoBit
.STAR_get_FASTA <- function(reference_path) {
    genome.fa <- file.path(reference_path, "resource", "genome.fa")
    if (!file.exists(genome.fa)) {
        genome.fa <- paste0(genome.fa, ".temp")
        genome.2bit <- file.path(reference_path, "resource", "genome.2bit")
        if (!file.exists(genome.2bit)) {
            .log(paste(genome.2bit, "not found"))
        }
        .log("Extracting temp genome FASTA from TwoBit file", "message")
        tmp <- rtracklayer::import(TwoBitFile(genome.2bit))
        rtracklayer::export(
            tmp,
            genome.fa, "fasta"
        )
        rm(tmp)
        gc()
    }
    return(genome.fa)
}

# Creates a temporary unzipped GTF for STAR
.STAR_get_GTF <- function(reference_path) {
    transcripts.gtf <- file.path(reference_path, "resource", "transcripts.gtf")
    if (!file.exists(transcripts.gtf)) {
        if (!file.exists(paste0(transcripts.gtf, ".gz"))) {
            .log(paste(paste0(transcripts.gtf, ".gz"), "not found"))
        }
        .log("Extracting temp Annotation GTF from GZ file", "message")
        R.utils::gunzip(paste0(transcripts.gtf, ".gz"), remove = FALSE,
            overwrite = TRUE)
        file.rename(transcripts.gtf, paste0(transcripts.gtf, ".temp"))
        transcripts.gtf <- paste0(transcripts.gtf, ".temp")
    }
    return(transcripts.gtf)
}

.STAR_clean_temp_FASTA_GTF <- function(reference_path) {
    .log("Cleaning temp genome / gene annotation files", "message")
    genome.fa <- file.path(reference_path, "resource", "genome.fa.temp")
    transcripts.gtf <- file.path(reference_path,
        "resource", "transcripts.gtf.temp")
    if (file.exists(genome.fa)) file.remove(genome.fa)
    if (file.exists(transcripts.gtf)) file.remove(transcripts.gtf)
}

.validate_STAR_fastq_samples <- function(fastq_1, fastq_2) {
    if (!is_valid(fastq_2)) {
        # assume single
        if (!all(file.exists(fastq_1))) {
            .log("Some fastq files were not found")
        }
    } else {
        if (length(fastq_2) != length(fastq_1)) {
            .log(paste("There must be equal numbers of",
                "forward and reverse fastq samples"))
        }
        if (!all(file.exists(fastq_1)) || !all(file.exists(fastq_2))) {
            .log("Some fastq files were not found")
        }
    }
    paired <- (length(fastq_1) == length(fastq_2))
    gzipped <- all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if (!gzipped &&
        (
            any(grepl(paste0("\\", ".gz", "$"), fastq_1)) ||
            (paired && any(grepl(paste0("\\", ".gz", "$"), fastq_2)))
        )
    ) {
        .log(paste("A mixture of gzipped and uncompressed",
            "fastq files found.", "You must supply either all",
            "gzipped or all uncompressed fastq files"))
    }
}

.validate_STAR_trim_sequence <- function(sequence) {
    if (length(sequence) != 1) {
        .log("Multiple adaptor sequences are not supported")
    }
    tryCatch({
        ACGT_sum <- sum(Biostrings::letterFrequency(
            Biostrings::DNAString(sequence),
            letters = "AGCT", OR = 0))
    }, error = function(e) ACGT_sum <- 0)
    if (nchar(sequence) != ACGT_sum) {
        .log("Adaptor sequence can only contain A, C, G or T")
    }
}