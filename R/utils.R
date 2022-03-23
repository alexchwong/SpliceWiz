# Exported miscellaneous functions

#' Converts genomic coordinates into a GRanges object
#'
#' This function takes a string vector of genomic coordinates and converts it
#' into a GRanges object.
#'
#' @details
#' Genomic coordinates can take one of the following syntax:
#' * `seqnames:start`
#' * `seqnames:start-end`
#' * `seqnames:start-end/strand`
#'
#' The following examples are considered valid genomic coordinates:
#' * "chr1:21535"
#' * "chr3:10550-10730"
#' * "X:51231-51330/-"
#' * "chrM:2134-5232/+"
#' @param coordinates A string vector of one or more genomic coordinates
#'   to be converted
#' @return A GRanges object that corresponds to the given coordinates
#' @examples
#' se <- SpliceWiz_example_NxtSE()
#'
#' coordinates <- rowData(se)$EventRegion
#'
#' gr <- CoordToGR(coordinates)
#' @md
#' @export
CoordToGR <- function(coordinates) {
    stopmsg <- paste(
        "Coordinates must take the form chrN:X, chrN:X-Y,",
        "chrN:X-Y/+ or chrN:X-Y/-"
    )
    coord_valid1 <- grepl(":", coordinates, fixed = TRUE)
    if(!all(coord_valid1)) {
        coord_error = coordinates[!coord_valid1]
        .log(paste(
            stopmsg,
            coord_error[1]
        ))
    }
    
    coord_hasdash <- grepl("-", 
        tstrsplit(coordinates, split = "/")[[1]], fixed = TRUE)
    coord_hasstrand <- grepl("/+", coordinates, fixed = TRUE) | 
        grepl("/-", coordinates, fixed = TRUE)

    if(any(!coord_hasdash)) {
        if(any(!coord_hasdash & !coord_hasstrand)) {
            coordinates[!coord_hasdash & !coord_hasstrand] <- 
                .coord_single_range_expand_nodash(
                    coordinates[!coord_hasdash & !coord_hasstrand])            
        }
        if(any(!coord_hasdash & coord_hasstrand)) {
            coordinates[!coord_hasdash & coord_hasstrand] <- 
                .coord_single_range_expand_withdash(
                    coordinates[!coord_hasdash & coord_hasstrand])      
        }
    }

    coord_hasstrand <- grepl("/+", coordinates, fixed = TRUE) | 
        grepl("/-", coordinates, fixed = TRUE)
        
    if(any(!coord_hasstrand)) {
        coordinates[!coord_hasstrand] <- paste0(coordinates[!coord_hasstrand], "/*")
    }
    
    # All coordinates should now have the correct format
    seqnames <- tstrsplit(coordinates, split = ":")[[1]]
    strands <- tstrsplit(coordinates, split = "/")[[2]]
    ranges <- tstrsplit(
        tstrsplit(
            coordinates, split = "/"
        )[[1]],
        split = ":"
    )[[2]]
    starts = as.numeric(tstrsplit(ranges, split = "-")[[1]])
    ends = as.numeric(tstrsplit(ranges, split = "-")[[2]])
    
    if(any(is.na(starts))) .log("Some coordinates have non-numeric starts")
    if(any(is.na(ends))) .log("Some coordinates have non-numeric ends")

    if(any(starts > ends)) .log("Some coordinates have negative ranges")

    return(GRanges(
        seqnames = seqnames, 
        ranges = IRanges(start = starts, end = ends),
        strand = strands
    ))
}

.coord_single_range_expand_withdash <- function(coords) {
    single_range <- as.numeric(
        tstrsplit(
            tstrsplit(coords, split = "/")[[1]],
            split = ":"
        )[[2]]
    )
    if(any(is.na(single_range))) .log(paste(
        "Some coordinates have single-nucleotide non-numeric ranges"))
    return(paste0(
        tstrsplit(
            tstrsplit(coords, split = "/")[[1]],
            split = ":"
        )[[1]], ":",
        paste(single_range, single_range, sep = "-"),
        "/", tstrsplit(coords, split = "/")[[2]]
    ))
}

.coord_single_range_expand_nodash <- function(coords) {
    single_range <- as.numeric(
        tstrsplit(coords, split = ":")[[2]]
    )
    if(any(is.na(single_range))) .log(paste(
        "Some coordinates have single-nucleotide non-numeric ranges"))
    return(paste0(
        tstrsplit(coords, split = ":")[[1]], ":",
        paste(single_range, single_range, sep = "-"),
        "/*"
    ))
}


#' Validates the given file as a valid COV file
#'
#' This function takes the path of a possible COV file and checks whether its
#' format complies with that of the COV format defined by this package.
#'
#' @details
#' COV files are BGZF-compressed files. The first 4 bytes of the file must
#' always be 'COV\1', distinguishing it from BAM or other files in BGZF format.
#' This function checks whether the given file complies with this.
#'
#' @param coverage_files A vector containing the file names of files to be
#'   checked
#' @return `TRUE` if all files are valid COV files. `FALSE` otherwise
#' @examples
#' se <- SpliceWiz_example_NxtSE()
#'
#' cov_files <- covfile(se)
#'
#' IsCOV(cov_files) # returns true if these are true COV files
#' @seealso [processBAM] [collateData]
#' @md
#' @export
IsCOV <- function(coverage_files) {
    for (i in coverage_files) {
        if (file.exists(i) && c_Check_Cov(normalizePath(i))) {
            # do nothing
        } else {
            return(FALSE)
        }
    }
    return(TRUE)
}
