#' SpliceWiz Example BAMs and NxtSE Experiment Object
#'
#' `SpliceWiz_example_bams()` is a wrapper function to obtain and make a local 
#' copy #' of 6 example files provided by ExperimentHub to
#' demonstrate the use of SpliceWiz. See [Example_data] for
#' a description of the provided BAM files. \cr\cr
#' `SpliceWiz_example_NxtSE()` retrieves a ready-made functioning
#' \linkS4class{NxtSE} object. The steps to reproduce this object is shown
#' in the example code in [MakeSE]
#'
#' @return
#' In `SpliceWiz_example_bams()`: returns a 2-column data frame containing
#'   sample names and BAM paths of the example dataset.
#'
#' In `SpliceWiz_example_NxtSE()`: returns a \linkS4class{NxtSE} object.
#' @examples
#'
#' # returns a data frame with the first column as sample names, and the
#' # second column as BAM paths
#'
#' SpliceWiz_example_bams()
#'
#' # Returns a NxtSE object created by the example bams aligned to the
#' # mock NxtSE reference
#'
#' se <- SpliceWiz_example_NxtSE()
#' @references
#' Generation of the mappability files was performed using SpliceWiz using
#' a method analogous to that described in:
#'
#' Middleton R, Gao D, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B,
#' Eyras E, Rasko JE, Ritchie W.
#' IRFinder: assessing the impact of intron retention on mammalian gene
#' expression.
#' Genome Biol. 2017 Mar 15;18(1):51.
#' \doi{10.1186/s13059-017-1184-4}
#' @name example-SpliceWiz-data
#' @aliases
#' SpliceWiz_example_bams SpliceWiz_example_NxtSE
#' @keywords package
#' @seealso [MakeSE]
#' @md
NULL

#' @describeIn example-SpliceWiz-data Returns a 2-column data frame, containing
#'   sample names and sample paths (in tempdir()) of example BAM files
#' @export
SpliceWiz_example_bams <- function() {
    bams <- example_bams()
    if (is.null(bams) || length(bams) != 6) stop("Example bam fetching failed")
    return(Find_Bams(tempdir()))
}

#' @describeIn example-SpliceWiz-data Returns a (in-memory / realized) NxtSE 
#' object that was pre-generated using the SpliceWiz example reference and 
#' example BAM files
#' @export
SpliceWiz_example_NxtSE <- function() {
    se <- readRDS(system.file("extdata",
        "example_NxtSE.Rds", package = "SpliceWiz"))
    covs <- Find_Samples(system.file("extdata", package = "SpliceWiz"), ".cov")
    covfile(se) <- covs$path
    se
}
