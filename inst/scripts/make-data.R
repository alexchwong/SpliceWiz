#' Sample NxtSE object using the example reference and example bam files
#'
#' This function reproduces the data stored in the inst/extdata directory.
#' This data includes:
#' 
#' * COV files produced by processing the example BAM dataset
#' * The NxtSE object produced by analysing the BAM dataset using the SpliceWiz
#'   workflow
#'
#' Note the Non-PolyA reference files (contained within 
#' "inst/extra-input-files") are from 
#' [IRFinder](https://github.com/williamritchie/IRFinder)
#' @examples
#' se <- readRDS(
#'   system.file("extdata", "example_NxtSE.Rds", package = "SpliceWiz")
#' )
#' se2 <- SpliceWiz_example_NxtSE()
#' identical(se, se2) # should return TRUE
#' @seealso 
#' [makeSE] for an outline of how to reproduce this data using the SpliceWiz
#' workflow
make_example_NxtSE <- function() {
    require(SpliceWiz)
    bams <- SpliceWiz_example_bams()
    buildRef(
        fasta = chrZ_genome(), gtf = chrZ_gtf(),
        reference_path = file.path(tempdir(), "Reference")
    )
    processBAM(
        bams$path, bams$sample,
        reference_path = file.path(tempdir(), "Reference"),
        output_path = file.path(tempdir(), "SpliceWiz_Output"),
        overwrite = TRUE, n_threads = 1
    )
    
    expr <- findSpliceWizOutput(file.path(tempdir(), "SpliceWiz_Output"))
    collateData(
        expr, 
        reference_path = file.path(tempdir(), "Reference"),
        output_path = file.path(tempdir(), "Collated_output")
    )
    
    se <- makeSE(collate_path = file.path(tempdir(), "Collated_output"))
    
    # Save COV files
    file.copy(
        expr$cov_file,
        file.path("../extdata", basename(expr$cov_file)),
        overwrite = TRUE
    )
    
    # De-identify COV files for validity once tempdir() no longer exists:
    covfile(se) <- rep("", 6)
    
    # Convert from HDF5-linked se to in-memory se:
    se <- realize_NxtSE(se, includeJunctions = TRUE)
    
    saveRDS(se, "../extdata/example_NxtSE.Rds")

    # Make example NxtSE with novel splice detection
    collateData(
        expr, 
        reference_path = file.path(tempdir(), "Reference"),
        output_path = file.path(tempdir(), "Collated_output_novel"),
        novelSplicing = TRUE
    )
    
    se <- makeSE(collate_path = file.path(tempdir(), "Collated_output_novel"))

    # De-identify COV files for validity:
    covfile(se) <- rep("", 6)
    
    # Convert from HDF5-linked se to in-memory se:
    se <- realize_NxtSE(se, includeJunctions = TRUE)
    
    saveRDS(se, "../extdata/example_NxtSE_novel.Rds")
}
