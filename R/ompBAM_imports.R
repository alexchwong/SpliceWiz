#' @useDynLib SpliceWiz, .registration = TRUE
#' @import Rcpp
#' @import zlibbioc
NULL

#' @export
idxstats <- function(bam_file, n_threads) {
    idxstats_pbam(bam_file, n_threads)
}
