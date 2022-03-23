#' Making Mappability Exclusion Region BED files from FASTA file
#'
#' Making a Mappability Region Exclusion BED file requires 3 steps.
#'   1) Generate mappability reads using `run_IRFinder_GenerateMapReads()`, 
#'      using the primary assembly genome FASTA as `genome.fa`
#'   2) Align `out.fa` to the corresponding genome using a genome splice-aware aligner such as STAR
#'   3) Process the aligned BAM file using `run_IRFinder_MapExclusionRegions()`
#' @examples
#' \dontrun{
#' FASTA = "/path/to/genome.fa"
#' run_IRFinder_GenerateMapReads(
#'      genome.fa = FASTA,
#'      out.fa = "/path/to/mappability_reads.fa",
#'      read_len = 70,
#'      read_stride = 10,
#'      error_pos = 35
#'  )
#'
#' # Now run STAR in the command line, e.g.:
#' # STAR \
#' # --genomeDir /path/to/Reference \
#' # --genomeLoad NoSharedMemory \
#' # --runThreadN 4 --outStd SAM --outSAMmode NoQS \
#' # --outSAMattributes None \
#' # --outFilterMultimapNmax 1 \
#' # --readFilesIn /path/to/mappability_reads.fa \
#' # > /path/to/genome_fragments.sam
#'
#' run_IRFinder_MapExclusionRegions(
#'      bamfile = "/path/to/genome_fragments.sam",
#'      output_file = "/path/to/MappabilityExclusionBED.txt",
#'      threshold = 4,
#'      includeCov = FALSE
#'  )
#' }
NULL

#' Sample NxtSE object using the NxtIRF example reference and example bam files
#'
#' This function generates an example reference, and runs IRFinder on example
#' bam files. This object is used in downstream examples throughout
#' the documentation of all functions that use NxtSE objects as input
#' @examples
#' se = readRDS(
#'   system.file("extdata", "example_NxtSE.Rds", package = "SpliceWiz")
#' )
#' se = SpliceWiz_example_NxtSE()
#' @seealso [makeSE()]
make_example_NxtSE <- function() {
    require(SpliceWiz)
    bams = SpliceWiz_example_bams()
    buildRef(
        fasta = chrZ_genome(), gtf = chrZ_gtf(),
        reference_path = file.path(tempdir(), "Reference")
    )
    processBAM(bams$path, bams$sample,
        reference_path = file.path(tempdir(), "Reference"),
        output_path = file.path(tempdir(), "SpliceWiz_Output"),
        overwrite = TRUE, n_threads = 1
    )
    expr = findSpliceWizOutput(file.path(tempdir(), "SpliceWiz_Output"))
    collateData(expr, 
      reference_path = file.path(tempdir(), "Reference"),
      output_path = file.path(tempdir(), "NxtIRF_output")
    )
    se = makeSE(collate_path = file.path(tempdir(), "NxtIRF_output"))
    
    # Save COV files
    file.copy(
        expr$cov_file,
        file.path("../extdata", basename(expr$cov_file)),
        overwrite = TRUE
    )
    
    # De-identify COV files for validity:
    covfile(se) <- rep("", 6)
    
    # Convert from HDF5-linked se to in-memory se:
    se <- realize_NxtSE(se)
    
    saveRDS(se, "../extdata/example_NxtSE.Rds")
}
