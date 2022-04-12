# Miscellaneous internal wrappers to SpliceWiz.dll

# Simple unzip function
# To check gunzip produces same output using SpliceWiz vs other utilities
run_Gunzip <- function(infile = "", outfile) {
    file_to_read <- normalizePath(infile)
    if (!file.exists(file_to_read)) {
        .log(paste("In run_Gunzip(),",
            file_to_read, "does not exist"))
    }
    if (!dir.exists(dirname(outfile))) {
        .log(paste("In run_Gunzip(),",
            dirname(outfile), "does not exist"))
    }
    c_gunzip(file_to_read, outfile)
}

# Gets a specific data frame in a gzipped multi-tabular text file
# If getting a small data frame situated at the beginning of a large file
#   e.g. in processBAM output, this is typically faster than data.table::fread
get_multi_DT_from_gz <- function(infile = "",
        block_headers = c("Header1", "Header2")) {
    file_to_read <- normalizePath(infile)
    if (!file.exists(file_to_read)) {
        .log(paste("In get_multi_DT_from_gz(),",
            file_to_read, "does not exist"))
    }
    df.list <- c_gunzip_DF(file_to_read, block_headers)
    for (i in seq_len(length(df.list))) {
        for (j in seq_len(length(df.list[[i]]))) {
            # suppressWarnings to supress "NAs introduced by coercion"
            # Intent of the function is to convert these to `NA`
            suppressWarnings({
                if (all(df.list[[i]][[j]] == "NA" |
                        !is.na(as.numeric(df.list[[i]][[j]])))) {
                    df.list[[i]][[j]] <- as.numeric(df.list[[i]][[j]])
                }
            })
        }
        df.list[[i]] <- as.data.table(df.list[[i]])
    }
    return(df.list)
}
