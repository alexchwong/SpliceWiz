#' Using Generalised linear models to analyse differential ASEs using edgeR
#'
#' These functions allow users to fit included/excluded counts using edgeR's
#' quasi-likelihood tests for differential Alternative Splice Events (ASEs)
#'
#' @details
#'
#' **edgeR** accounts appropriately for zero-counts which are often problematic
#' as PSI approaches zero or one, leading to false positives.
#' The following functions allow users to define model formulas to test relative
#' expressions of included / excluded counts (to assess whether isoforms are
#' differentially regulated, in isolation), as well as together as an 
#' interaction (the latter provides results of differential ASE analysis)
#'
#' See the examples section for a brief explanation of how to use these 
#' functions.
#'
#' See also [ASE-methods] for further explanations of results output.
#'
#' @param se The \linkS4class{NxtSE} object created by `makeSE()`. To reduce
#'   runtime and avoid excessive multiple testing, consider filtering
#'   the object using [applyFilters]
#' @param strModelFormula A string specifying the model formula to fit isoform
#'   counts to assess differential expression in isolation. Should take the form
#'   of `"~0 + batch1 + batch2 + test_factor"`, where `batch1` and `batch2` are
#'   batch factors (if any), and `test_factor` is the variate of interest.
#' @param strASEFormula A string specifying the model formula to fit PSIs 
#'   (isoform ratios). The variate of interest should be specified as an
#'   interactiion term with `ASE`. For example, following the above example, 
#'   the ASE formula should be 
#'   `"~0 + batch1 + batch2 + test_factor + test_factor:ASE"`
#' @param model_IncExc A model matrix in which to model differential expression
#'   of isoform counts in isolation. The number of rows must equal that of the
#'   number of samples in `se`
#' @param model_ASE A model matrix in which to model differential PSIs.
#'   The number of rows must be twice that of the number of samples in `se`, the
#'   first half are for included counts, and the second half are for excluded
#'   counts. See example below.
#' @param IRmode (default `all`) Choose the approach to quantify IR events.
#'   Default `all` considers all introns as potentially retained, and calculates
#'   IR-ratio based on total splicing across the intron using the "SpliceOver"
#'   or "SpliceMax" approach (see [collateData]). Other options include 
#'   `annotated` which calculates IR-ratios for annotated introns only, and
#'   `annotated_binary` which calculates PSI considering the "included"
#'   isoform as the IR-transcript, and the "excluded" transcript is
#'   quantified from splice counts only across the exact intron 
#'   (but not that of overlapping introns). IR-ratio are denoted as "IR" events,
#'   whereas PSIs calculated using IR and intron-spliced binary alternatives are
#'   denoted as "RI" events.
#' @param filter_antiover,filter_antinear Whether to remove novel IR events that
#'   overlap over or near anti-sense genes. Default will exclude antiover but
#'   not antinear introns. These are ignored if strand-specific RNA-seq 
#'   protocols are used.
#' @param conditionList A list (or vector) of condition values of which to
#' calculate mean PSIs
#' @param condition The name of the column containing the condition values in
#'   `colData(se)`
#' @param fit The output returned by the `fitASE_edgeR` and 
#'   `fitASE_edgeR_custom` functions.
#' @param coef_IncExc,coef_ASE model coefficients to be dropped for LRT test
#'   between full and reduced models. Directly parsed onto `edgeR::glmQLFTest`.
#'   See `?edgeR::glmQLFTest` for details
#' @param contrast_IncExc,contrast_ASE numeric vector specifying one or more #'   contrasts of the linear model coefficients to be tested.
#'   Directly parsed onto `edgeR::glmQLFTest`. See `?edgeR::glmQLFTest` for
#'   details
#' @param results The return value of `testASE_edgeR()`, to be used as input to
#'   append mean and delta PSI values onto.
#' @param useQL (default `TRUE`) Whether to use edgeR's quasi-likelihood method
#'   to help reduce false positives from near-zero junction / intron counts.
#'   NB: edgeR's quasi-likelihood method is run with legacy method 
#'   (Lun and Smyth (2017)).
#' @return 
#'   `fitASE_edgeR` and `fitASE_edgeR_custom` returns a named list containing
#'   the following:
#'   * `IncExc`, `ASE`: `DGEGLM` objects containing the fitted models for 
#'     isoform counts and PSIs, respectively
#'   * `model_IncExc`, `model_ASE`: model matrices of the above fitted models.
#' 
#'   `testASE_edgeR()` returns a data.table containing the following:
#'   * EventName: The name of the ASE event. This identifies each ASE
#'     in downstream functions including [makeMeanPSI], [makeMatrix],
#'     and [plotCoverage]
#'   * EventType: The type of event. See details section above.
#'   * EventRegion: The genomic coordinates the event occupies. This spans the
#'     most upstream and most downstream splice junction involved in the ASE,
#'     and is use to guide the [plotCoverage] function.
#'   * flags: Indicates which isoforms are NMD substrates and/or which are
#'     formed by novel splicing only.
#'
#'   **edgeR specific output** equivalent to statistics returned by 
#'   `edgeR::topTags()`:
#'   * logFC, logCPM, F, PValue, FDR: log fold change, log counts per million,
#'     F statistic, p value and (Benjamini Hochberg) adjusted p values of the
#'     differential PSIs for the contrasts or coefficients tested.
#'   * inc/exc_(...): edgeR statistics corresponding to 
#'     differential expression testing for raw included / excluded counts
#'     in isolation (not of the PSIs).
#'
#'   `addPSI_edgeR()` appends the following columns to the above output
#'   * AvgPSI_X: the average percent spliced in / percent
#'     IR levels for condition X. Note this is a
#'     geometric mean, based on the arithmetic mean of logit PSI values.
#'   * deltaPSI: The difference in PSI between the mean values of the two
#'     conditions.
#'   * abs_deltaPSI: The absolute value of difference in PSI between 
#'     the mean values of the two conditions.
#'
#' @examples
#' # Load the NxtSE object and set up the annotations
#' # - see ?makeSE on example code of generating this NxtSE object
#' se <- SpliceWiz_example_NxtSE()
#'
#' colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' require("edgeR")
#'
#' fit <- fitASE_edgeR(
#'     se, 
#'     strModelFormula = "~0 + replicate + treatment", 
#'     strASEFormula = "~0 + replicate + treatment + treatment:ASE"
#' )
#' 
#' # Get coefficient terms of Included / Excluded counts isolated model
#' colnames(fit$model_IncExc)
#' # [1] "replicateP" "replicateQ" "replicateR" "treatmentB"
#' 
#' # Get coefficient terms of PSI model
#' colnames(fit$model_ASE)
#' # [1] "replicateP" "replicateQ" "replicateR" "treatmentB"            
#' # [5] "treatmentA:ASEIncluded" "treatmentB:ASEIncluded"
#' 
#' # Contrast between treatment "B" against treatment "A"
#' res <- testASE_edgeR(se, fit,
#'     contrast_IncExc = c(0,0,0,1),
#'     contrast_ASE = c(0,0,0,0,-1,1)
#' )
#'
#' ### # Add mean PSI values to results:
#' res_withPSI <- addPSI_edgeR(res, se, "treatment", c("B", "A"))
#'
#'
#' ### Using custom model matrices to model counts
#' #   - the equivalent analysis can be performed as follows:
#'
#' # Sample annotations for isoform count expressions
#' colData <- as.data.frame(colData(se))
#'
#' # Sample annotations for isoform count PSI analysis
#' colData_ASE <- rbind(colData, colData)
#' colData_ASE$ASE <- rep(c("Included", "Excluded"), each = nrow(colData))
#' rownames(colData_ASE) <- c(
#'     paste0(rownames(colData), ".Included"),
#'     paste0(rownames(colData), ".Excluded")
#' )
#'
#' model_IncExc <- model.matrix(
#'     ~0 + replicate + treatment,
#'     data = colData
#' )
#'
#' model_ASE <- model.matrix(
#'     ~0 + replicate + treatment + treatment:ASE,
#'     data = colData_ASE
#' )
#' 
#' fit <- fitASE_edgeR_custom(se, model_IncExc, model_ASE)
#'
#' res_customModel <- testASE_edgeR(se, fit,
#'     contrast_IncExc = c(0,0,0,1),
#'     contrast_ASE = c(0,0,0,0,-1,1)
#' )
#'
#' # Check this produces identical results:
#' identical(res_customModel, res)
#'
#' ### Time series examples using edgeR and splines 
#' # - similar to section 4.8 in the edgeR vignette
#'
#' colData(se)$timepoint <- rep(c(1,2,3), each = 2)
#' colData(se)$batch <- rep(c("1", "2"), 3)
#'
#' # First, we set up a polynomial spline with 2 degrees of freedom:
#' Time <- poly(colData(se)$timepoint, df = 2)
#'
#' # Next, we define the batch factor:
#' Batch <- factor(colData(se)$batch)
#'
#' # Finally, we construct the same factors for ASE analysis. Note that
#' #   each factor must be repeated twice
#'
#' Time_ASE <- rbind(Time, Time)
#' Batch_ASE <- c(Batch, Batch)
#' ASE <- factor(
#'     rep(c("Included", "Excluded"), each = nrow(colData(se)))
#' )
#'
#' # Now, we set up the model matrices for isoform and PSI count modelling
#' model_IncExc <- model.matrix(~0 + Batch + Time)
#' model_ASE <- model.matrix(~0 + Batch_ASE + Time_ASE + Time_ASE:ASE)
#'
#' fit <- fitASE_edgeR_custom(se, model_IncExc, model_ASE)
#'
#' # Note the coefficients of interest in the constructed models:
#'
#' colnames(model_IncExc)
#' # [1] "Batch1" "Batch2" "Time1"  "Time2" 
#'
#' colnames(model_ASE)
#' # [1] "Batch_ASE1" "Batch_ASE2" "Time_ASE1" "Time_ASE2"
#' # [5] "Time_ASE1:ASEIncluded" "Time_ASE2:ASEIncluded"
#' 
#' # We are interested in a model in which `Time` is excluded, thus:
#' 
#' res <- testASE_edgeR(se, fit,
#'     coef_IncExc = 3:4,
#'     coef_ASE = 5:6
#' )
#'
#' # Finally, add PSI values for each time point:
#'
#' res_withPSI <- addPSI_edgeR(res, se, "timepoint", c(1, 2, 3))
#'
#' @name ASE-GLM-edgeR
#' @references
#' Lun A, Smyth G (2017).
#' 'No counts, no variance: allowing for loss of degrees of freedom when
#' assessing biological variability from RNA-seq data' 
#' Stat Appl Genet Mol Biol, 017 Apr 25;16(2):83-93.
#' \url{https://doi.org/10.1515/sagmb-2017-0010}
NULL

#' @describeIn ASE-GLM-edgeR Use edgeR to fit counts and ASE models with a
#'   given design formula
#' @export
fitASE_edgeR <- function(
    se,
    strModelFormula, strASEFormula,
    useQL = TRUE,
    IRmode = c("all", "annotated", "annotated_binary"),
    filter_antiover = TRUE, filter_antinear = FALSE
) {
    .check_package_installed("edgeR", "3.32.0")
    models <- .ASE_get_models(se, strModelFormula, strASEFormula)
    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)

    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")

    .log("Fitting edgeR contrast for included / excluded counts separately",
        "message")
    fit_IncExc <- .ASE_edgeR_fitIncExc(se_use, models$design)

    rownames_Inc <- rownames(fit_IncExc$fit)[
        grepl(".Included", rownames(fit_IncExc$fit))]
    rownames_Inc <- sub(".Included", "", rownames_Inc)
    rownames_Exc <- rownames(fit_IncExc$fit)[
        grepl(".Excluded", rownames(fit_IncExc$fit))]
    rownames_Exc <- sub(".Excluded", "", rownames_Inc)

    se_use <- se_use[rowData(se_use)$EventName %in% rownames_Inc &
        rowData(se_use)$EventName %in% rownames_Exc,]

    .log("Fitting edgeR contrast for included / excluded counts together",
        "message")
    fit_ASE <- .ASE_edgeR_fitASE(se_use, models$design_ASE, useQL)

    return(list(
        IncExc = fit_IncExc$fit,
        ASE = fit_ASE$fit,
        model_IncExc = fit_IncExc$model,
        model_ASE = fit_ASE$model
    ))
}

#' @describeIn ASE-GLM-edgeR Use edgeR to fit counts and ASE models with a
#'   given design formula
#' @export
fitASE_edgeR_custom <- function(
        se, model_IncExc, model_ASE,
        useQL = TRUE,
        IRmode = c("all", "annotated", "annotated_binary"),
        filter_antiover = TRUE, filter_antinear = FALSE
) {
    .check_package_installed("edgeR", "3.32.0")
    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)

    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")

    .log("Fitting edgeR contrast for included / excluded counts separately",
        "message")
    fit_IncExc <- .ASE_edgeR_fitIncExc(se_use, model_IncExc)

    rownames_Inc <- rownames(fit_IncExc$fit)[
        grepl(".Included", rownames(fit_IncExc$fit))]
    rownames_Inc <- sub(".Included", "", rownames_Inc)
    rownames_Exc <- rownames(fit_IncExc$fit)[
        grepl(".Excluded", rownames(fit_IncExc$fit))]
    rownames_Exc <- sub(".Excluded", "", rownames_Inc)

    se_use <- se_use[rowData(se_use)$EventName %in% rownames_Inc &
        rowData(se_use)$EventName %in% rownames_Exc,]

    .log("Fitting edgeR contrast for included / excluded counts together",
        "message")
    fit_ASE <- .ASE_edgeR_fitASE(se_use, model_ASE, useQL)

    return(list(
        IncExc = fit_IncExc$fit,
        ASE = fit_ASE$fit,
        model_IncExc = fit_IncExc$model,
        model_ASE = fit_ASE$model
    ))
}

#' @describeIn ASE-GLM-edgeR Use edgeR to return differential ASE results. coef
#'   and contrast are parsed onto edgeR's glmQLFTest function
#' @export
testASE_edgeR <- function(
    se, fit,
    coef_IncExc = ncol(fit[["model_IncExc"]]), 
    contrast_IncExc = NULL,
    coef_ASE = ncol(fit[["model_ASE"]]), 
    contrast_ASE = NULL
) {
    useQL <- !is.null(fit$ASE$df.residual.zeros)
    
    se_use <- se[rownames(se) %in% rownames(fit$ASE)]
    
    if(useQL) {
        qlf_IncExc <- edgeR::glmQLFTest(
            fit$IncExc,
            coef = coef_IncExc, contrast = contrast_IncExc
        )    
    } else {
        qlf_IncExc <- edgeR::glmLRT(
            fit$IncExc,
            coef = coef_IncExc, contrast = contrast_IncExc
        )  
    }
    
    res_IncExc <- edgeR::topTags(qlf_IncExc, n = nrow(fit$fit_IncExc))
    res_IncExc$table$EventName <- rownames(res_IncExc$table)

    res_IncExc <- as.data.table(res_IncExc$table)
    res.inc <- res_IncExc[grepl(".Included", get("EventName"))]
    res.inc[, c("EventName") :=
        sub(".Included","",get("EventName"), fixed=TRUE)]
    res.exc <- res_IncExc[grepl(".Excluded", get("EventName"))]
    res.exc[, c("EventName") :=
        sub(".Excluded","",get("EventName"), fixed=TRUE)]
        
    if(useQL) {
        qlf_ASE <- edgeR::glmQLFTest(
            fit$ASE,
            coef = coef_ASE, contrast = contrast_ASE
        )
    } else {
        qlf_ASE <- edgeR::glmLRT(
            fit$ASE,
            coef = coef_ASE, contrast = contrast_ASE
        )
    }
    res_ASE <- edgeR::topTags(qlf_ASE, n = nrow(fit$fit_IncExc))
    res_ASE$table$EventName <- rownames(res_ASE)
    
    res_ASE <- as.data.table(res_ASE$table)
    colnames(res.inc)[-ncol(res.inc)] <- paste(
        "Inc", colnames(res.inc)[-ncol(res.inc)], sep = ".")
    colnames(res.exc)[-ncol(res.exc)] <- paste(
        "Exc", colnames(res.exc)[-ncol(res.exc)], sep = ".")
    
    res_ASE <- res_ASE[res.inc, on = "EventName"]
    res_ASE <- res_ASE[res.exc, on = "EventName"]
    
    orderCol <- ifelse(useQL, "F", "LR")
    setorderv(res_ASE, orderCol, order = -1)
    
    rowData <- as.data.frame(rowData(se_use))
    rowData.DT <- .ASE_add_flags(as.data.table(rowData[,
        c("EventName","EventType","EventRegion", "NMD_direction")]))
    res <- cbind(
        res_ASE[,c("EventName")],
        res_ASE[,-c("EventName")]
    )
    res <- rowData.DT[res, on = "EventName"]
    return(res)
}

#' @describeIn ASE-GLM-edgeR Adds average and delta PSIs of conditions of
#'   interest onto results produced by testASE_edgeR(). Note this is done
#'   automatically for other methods described in `ASE-methods`.
#' @export
addPSI_edgeR <- function(
    results, se,
    condition, conditionList
) {
    if(!all(c("EventName","EventType","EventRegion", "flags") %in%
        colnames(results)
    )) {
        .log("`results` must be the output of testASE_edgeR() function")
    }
    if(!(condition %in% colnames(colData(se)))) {
        .log("`condition` must be a column name in colData(se)")
    }
    conds <- unique(colData(se)[, condition])
    if(!all(conditionList %in% conds)) {
        .log(paste("All elements `conditionList` must be elements of",
            "colData(se)[, condition]"
        ))
    }
    results <- .ASE_add_diag_multi(results, se, condition, conditionList)
    return(results)
}

.ASE_get_models <- function(se, 
        strModelFormula, strASEFormula
) {
    colData <- as.data.frame(colData(se))
    design <- model.matrix(
        as.formula(strModelFormula),
        data = colData
    )
    colData_ASE <- rbind(colData, colData)
    colData_ASE$ASE <- rep(c("Included", "Excluded"), each = nrow(colData))
    rownames(colData_ASE) <- c(
        paste0(rownames(colData), ".Included"),
        paste0(rownames(colData), ".Excluded")
    )
    design_ASE <- model.matrix(
        as.formula(strASEFormula),
        data = colData_ASE
    )
    return(list(
        design = design,
        design_ASE = design_ASE,
        colData = colData,
        colData_ASE = colData_ASE
    ))
}

.ASE_edgeR_fitIncExc <- function(se, model) {
    countData <- as.matrix(rbind(assay(se, "Included"),
        assay(se, "Excluded")))
    rowData <- as.data.frame(rowData(se))
    colData <- colData(se)
    rownames(colData) <- colnames(se)
    colnames(countData) <- rownames(colData)
    rownames(countData) <- c(
        paste(rowData$EventName, "Included", sep="."),
        paste(rowData$EventName, "Excluded", sep=".")
    )

    y <- edgeR::DGEList(counts=countData, remove.zeros = FALSE)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, model)
    
    fit <- edgeR::glmQLFit(y, model, legacy = TRUE)
    return(list(
        fit = fit,
        model = model
    ))
}

.ASE_edgeR_fitASE <- function(se, model, useQL = TRUE) {
    countData <- as.matrix(cbind(assay(se, "Included"),
        assay(se, "Excluded")))

    rowData <- as.data.frame(rowData(se))
    colData <- as.data.frame(colData(se))
    colData <- rbind(colData, colData)
    rownames(colData) <- c(
        paste(colnames(se), "Included", sep="."),
        paste(colnames(se), "Excluded", sep=".")
    )
    colData$ASE <- rep(c("Included", "Excluded"), each = ncol(se))
    colnames(countData) <- rownames(colData)
    rownames(countData) <- rowData$EventName

    y <- edgeR::DGEList(counts=countData, remove.zeros = FALSE)
    y <- edgeR::estimateDisp(y, model)
    y$offset <- 1
    
    if(useQL) {
        fit <- edgeR::glmQLFit(y, model, legacy = TRUE)    
    } else {
        fit <- edgeR::glmFit(y, model)
    }
    return(list(
        fit = fit,
        model = model
    ))
}
