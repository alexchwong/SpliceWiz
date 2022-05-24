.onAttach <- function(libname, pkgname) {
    SW_threads <- setSWthreads(0)
    packageStartupMessage("SpliceWiz package loaded with ", SW_threads, 
        " threads")
    packageStartupMessage(
        "Use setSWthreads() to set the number of SpliceWiz threads")
    invisible()
}

#' Sets the number of threads used by SpliceWiz
#' 
#' SpliceWiz uses the computationally efficient packages `fst` and `data.table`
#' to compute file and data operations, respectively. Both packages make use
#' of parallelisation. If excessive number of threads are allocated, it may
#' impact the running of other operations on your system. Use this function to
#' manually allocate the desired number of threads
#' @param threads (default `0`) The number of threads for SpliceWiz to use.
#'   Set as `0` to use the recommended number of threads appropriate for the
#'   system (approximately half the available threads)
#' @return Nothing.
#' @examples
#' setSWthreads(0)
#' @export
setSWthreads <- function(threads = 0) {
    system_threads <- parallel::detectCores()
    DT_threads <- getDTthreads()
    FST_threads <- threads_fst()

    if(threads == 0) {
        if(DT_threads <= ceiling(system_threads / 2)) {
            SW_threads <- DT_threads
        } else {
            SW_threads <- ceiling(system_threads / 2)
        }
    } else if(threads <= system_threads) {
        SW_threads <- threads
    } else {
        .log("Requested threads exceed system resources")
    }
    
    setDTthreads(SW_threads)
    threads_fst(SW_threads)
    return(SW_threads)
}

.getSWthreads <- function() {
    DT_threads <- getDTthreads()
    FST_threads <- threads_fst()
    return(min(DT_threads, FST_threads))
}