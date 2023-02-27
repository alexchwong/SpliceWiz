#' SpliceWiz Example BAMs and NxtSE Experiment Object
#'
#' `SpliceWiz_example_bams()` is a wrapper function to obtain and make a local 
#' copy
#' of 6 example files provided by the NxtIRFdata companion package to
#' demonstrate the use of SpliceWiz. See [NxtIRFdata::example_bams] for
#' a description of the provided BAM files. \cr\cr
#' `SpliceWiz_example_NxtSE()` retrieves a ready-made functioning
#' \linkS4class{NxtSE} object. The steps to reproduce this object is shown
#' in the example code in [makeSE]
#'
#' @param novelSplicing Whether to import an example NxtSE with novel splice
#'   event discovery.
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
#' @seealso [makeSE]
#' @md
NULL

#' @describeIn example-SpliceWiz-data Returns a 2-column data frame, containing
#'   sample names and sample paths (in tempdir()) of example BAM files
#' @export
SpliceWiz_example_bams <- function() {
    bams <- NxtIRFdata::example_bams()
    if (is.null(bams) || length(bams) != 6) stop("Example bam fetching failed")
    return(
        data.frame(
            sample = tstrsplit(basename(bams), split = ".", fixed = TRUE)[[1]],
            path = bams
        )
    )
}

#' @describeIn example-SpliceWiz-data Returns a (in-memory / realized) NxtSE 
#' object that was pre-generated using the SpliceWiz example reference and 
#' example BAM files
#' @export
SpliceWiz_example_NxtSE <- function(novelSplicing = FALSE) {
    if(!novelSplicing) {
        se <- readRDS(system.file("extdata",
            "example_NxtSE.Rds", package = "SpliceWiz"))
    } else {
        se <- readRDS(system.file("extdata",
            "example_NxtSE_novel.Rds", package = "SpliceWiz"))    
    }
    covs <- findSamples(system.file("extdata", package = "SpliceWiz"), ".cov")
    se@metadata[["sourcePath"]] <- dirname(covs$path[1])
    covfile(se) <- covs$path
    se
}
