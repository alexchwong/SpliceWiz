.onAttach <- function(libname, pkgname) {

    # Detect optional dependencies, prompt to install if missing:
    CRAN_deps <- c("DoubleExpSeq", "egg", "DBI", "crayon")
    Bioc_deps <- c("DESeq2", "limma", "satuRn", "edgeR", 
        "GO.db", "fgsea", "Rsubread")

    CRAN_deps_inst <- vapply(CRAN_deps, .check_package_installed,
        FUN.VALUE = logical(1), "0.0.0", "silent")
    Bioc_deps_inst <- vapply(Bioc_deps, .check_package_installed,
        FUN.VALUE = logical(1), "0.0.0", "silent")
    
    if(!all(CRAN_deps_inst) | !all(Bioc_deps_inst)) {
        CRAN_deps_install <- c()
        Bioc_deps_install <- c()
            
        CRAN_deps_prompt <- CRAN_deps[!CRAN_deps_inst]
        if(length(CRAN_deps_prompt) > 0) {
            CRAN_deps_prompt <- paste0('"', CRAN_deps_prompt, '"')

            CRAN_deps_install <- paste0(
                "install.packages(", 
                ifelse(length(CRAN_deps_prompt) > 1, "c(", ""),
                paste(CRAN_deps_prompt, collapse = ", "),
                ifelse(length(CRAN_deps_prompt) > 1, "))", ")")
            )        
        }
        
        Bioc_deps_prompt <- Bioc_deps[!Bioc_deps_inst]
        if(length(Bioc_deps_prompt) > 0) {
            Bioc_deps_prompt <- paste0('"', Bioc_deps_prompt, '"')

            Bioc_deps_install <- paste0(
                "BiocManager::install(", 
                ifelse(length(Bioc_deps_prompt) > 1, "c(", ""),
                paste(Bioc_deps_prompt, collapse = ", "),
                ifelse(length(Bioc_deps_prompt) > 1, "))", ")")
            )
        }
        
        .log(paste0(
            "Some optional dependencies are not installed.\n",
            "For full SpliceWiz functionally, install all dependencies",
            " by running the following:\n\n",
            CRAN_deps_install, ifelse(is.null(CRAN_deps_install), "", "\n\n"),
            Bioc_deps_install, ifelse(is.null(Bioc_deps_install), "", "\n\n")
        ), "message")
    }
    
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