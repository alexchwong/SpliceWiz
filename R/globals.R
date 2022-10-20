# Global Internal functions

globalVariables(c(":=", "."))

collateData_version <- "0.99.5"

buildRef_version <- "0.99.5"

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

is_valid <- function(x) {
    !missing(x) &&
    !is.null(x) && 
    length(x) > 0 &&
    (
        length(x) > 1 || 
        isS4(x) || 
        is(x, "data.frame") || 
        !is.na(x)
    ) &&
    (
        length(x) > 1 || 
        !is.character(x) || 
        (x != "" && x != "(none)")
    )
}

.log <- function(msg = "",
        type = c("error", "warning", "silent", "message"),
        use_system_time = TRUE,
        ...
) {
    type <- match.arg(type)
    if (use_system_time) {
        msg <- paste(format(Sys.time(), "%b %d %X"), msg)
    }
    if (type == "error") {
        stop(msg, call. = FALSE)
    } else if (type == "warning") {
        warning(msg, call. = FALSE)
    } else {
        message(msg, ...)
    }
}

.check_package_installed <- function(
        package = "DESeq2", version = "0.0.0",
        returntype = c("error", "warning", "silent")
) {
    res <- tryCatch(
        ifelse(packageVersion(package) >= version, TRUE, FALSE),
        error = function(e) FALSE)
    if (!res) {
        returntype <- match.arg(returntype)
        .log(paste(package, "version", version, "is not installed;",
            "and is required for this function"), type = returntype)
    }
    return(res)
}

.colourise <- function(text, color) {
    if (!.check_package_installed("crayon", returntype = "silent") ||
        !crayon::has_color()
    ) return(text)
    return(crayon::style(text, color))
}

.nxtcat <- function(...) {
    cat(sprintf(...))
}

.validate_threads <- function(n_threads, as_BPPARAM = TRUE, 
        useSnowParam = FALSE,
        ...
) {
    n_threads_to_use <- as.numeric(n_threads)
    if (is.na(n_threads_to_use)) {
        .log("n_threads must be a numeric value")
    }
    if (n_threads_to_use > parallel::detectCores()) {
        n_threads_to_use <- max(1, parallel::detectCores())
    }
    if (as_BPPARAM) {
        if (useSnowParam || Sys.info()["sysname"] == "Windows") {
            BPPARAM_mod <- BiocParallel::SnowParam(n_threads_to_use, ...)
            .log(paste("Using SnowParam", BPPARAM_mod$workers, "threads"),
                "message")
            setSWthreads(1) 
            # SnowParam doesn't count as a fork for data.table or fst          
        } else {
            BPPARAM_mod <- BiocParallel::MulticoreParam(n_threads_to_use, ...)
            .log(paste("Using MulticoreParam", BPPARAM_mod$workers, "threads"),
                "message")
        }
        return(BPPARAM_mod)
    } else {
        return(n_threads_to_use)
    }
}

.restore_threads <- function(n_threads) {
    setSWthreads(n_threads)
}

.split_vector <- function(vector = "", n_workers = 1) {
    if (!is.numeric(n_workers) || n_workers < 1)
        .log("n_workers must be at least 1")

    n_workers_use <- as.integer(n_workers)
    if (length(vector) < 1) .log("vector to split must be of length at least 1")

    if (n_workers_use > length(vector)) n_workers_use <- length(vector)
    vector_starts <- round(seq(1, length(vector) + 1,
        length.out = n_workers_use + 1))
    vector_starts <- unique(vector_starts)

    return_val <- list()
    for (i in seq_len(length(vector_starts) - 1)) {
        return_val[[i]] <- vector[
            seq(vector_starts[i], vector_starts[i + 1] - 1)]
    }
    return(return_val)
}

# Semi-join a data.table. Equivalent to dplyr::semi_join(A, B, by = by)
.semi_join_DT <- function(A, B, by, nomatch = 0) {
    A[A[B, on = by, which = TRUE, nomatch = nomatch]]
}

# Converts data table to GRanges object, preserving info
.grDT <- function(DT, ...) {
    if (nrow(DT) == 0) return(GenomicRanges::GRanges(NULL))
    makeGRangesFromDataFrame(as.data.frame(DT), ...)
}

.grlGaps <- function(grl) {
    psetdiff(unlist(range(grl), use.names = TRUE), grl)
}

.make_path_relative <- function(files, relative_to) {
    if (is.null(files)) return(files)
    if (!all(file.exists(files))) .log("Some files do not exist")
    if (!all(file.exists(relative_to))) .log("Some directories do not exist")
    if (length(relative_to) == 1) relative_to <- rep(relative_to, length(files))

    if (Sys.info()["sysname"] == "Windows") {
        files <- normalizePath(files, winslash = "/")
        relative_to <- normalizePath(relative_to, winslash = "/")
    } else {
        files <- normalizePath(files)
        relative_to <- normalizePath(relative_to)
    }
    out <- c()
    for (i in seq_len(length(files))) {
        f <- files[i]
        base <- relative_to[i]

        common <- sub("^([^|]*)[^|]*(?:\\|\\1[^|]*)$", "^\\1/?",
            paste0(base, "|", f))

        out <- c(out, paste0(gsub("[^/]+/?", "../", sub(common, "", base)),
            sub(common, "", f)))
    }
    return(out)
}

# Compatibility for running inside a shiny withProgress block
dash_progress <- function(message = "", total_items = 1, add_msg = FALSE) {
    if (total_items != round(total_items) | total_items < 1) {
        .log("dash_progress needs at least 1 item")
    }
    if (add_msg) {
        .log(message, "message")
    }
    has_shiny <- .check_package_installed(
        package = "shiny", returntype = "silent")
    if (has_shiny) {
        session <- shiny::getDefaultReactiveDomain()
        if (!is.null(session)) {
            shiny::incProgress(1 / total_items, message = message)
        }
    }
}

# GUI specific functions

update_data_frame <- function(existing_df, new_df) {
    # add extra samples to existing df
    DT1 <- as.data.table(existing_df)
    DT2 <- as.data.table(new_df)

    common_cols <- intersect(names(DT1)[-1], names(DT2)[-1])
    new_cols <- names(DT2)[!(names(DT2) %in% names(DT1))]

    if(!all(DT2$sample %in% DT1$sample)) {
        DT_add <- DT2[!(sample %in% DT1$sample)]
        if(length(new_cols) > 0) DT_add <- DT_add[, c(new_cols) := NULL]
        newDT <- rbind(DT1, DT_add, fill = TRUE)
    } else {
        newDT <- copy(DT1)
    }

    if(length(new_cols) > 0) {
        DT_tomerge <- copy(DT2)
        if(length(common_cols) > 0) {
            DT_tomerge[, c(common_cols) := NULL]
        }
        newDT <- merge(newDT, DT_tomerge, all = TRUE, by = "sample")
    }

    # now update conflicting values
    if(length(common_cols) > 0 & any(DT2$sample %in% DT1$sample)) {
        DT_toupdate <- DT2[(sample %in% DT1$sample)]
        if(length(new_cols) > 0) {
            DT_toupdate <- DT_toupdate[, c(new_cols) := NULL]
        }
        newDT[DT_toupdate, on = .(sample), 
            (common_cols) := mget(paste0("i.", common_cols))]
    }
    return(as.data.frame(newDT))
}

