#' Differential Alternative Splicing Event analysis
#'
#' Use Limma, DESeq2, DoubleExpSeq and satuRn wrapper functions to test for
#' differential Alternative Splice Events (ASEs)
#'
#' @details
#'
#' Using **limma**, SpliceWiz models included and excluded counts as log-normal
#' distributed, whereas
#' using **DESeq2**, SpliceWiz models included and excluded counts as negative
#' binomial distributed with dispersion shrinkage according to their mean count
#' expressions.
#' For **limma** and **DESeq2**, differential ASE are considered as the
#' "interaction" between included and excluded splice counts for each sample.
#'  See [this vignette](https://rpubs.com/mikelove/ase) for an explanation of
#' how this is done.
#'
#' SpliceWiz's **limma** wrapper implements an additional filter where ASEs with
#' an average cpm values of either Included or Excluded counts are less than 1.
#' **DESeq2** has its own method for handling outliers, which seems to work well
#' for handling situations where PSI ~= 0 or PSI ~= 1.
#'
#' Time series are supported by SpliceWiz to a limited extent. Time series
#' analysis is only performed via DESeq2 (using its "LRT" mode). To activate
#' time series differential analysis, run `ASE_DESeq()` specifying `test_factor`
#' as the column of numeric values containing time series data. The `test_nom`
#' and `test_denom` parameters must be left blank. See example below.
#'
#' Using **DoubleExpSeq**, included and excluded counts are modeled using
#' the generalized beta prime distribution, using empirical Bayes shrinkage
#' to estimate dispersion.
#'
#' Using **satuRn**, included and excluded counts are modeled using
#' the quasi-binomial distribution in a generalised linear model.
#'
#' **EventType** are as follow:
#' * `IR` = intron retention (IR-ratio) - all introns are considered
#' * `MXE` = mutually exclusive exons
#' * `SE` = skipped exons
#' * `AFE` = alternate first exon
#' * `ALE` = alternate last exon
#' * `A5SS` = alternate 5'-splice site
#' * `A3SS` = alternate 3'-splice site
#' * `RI` = (known / annotated) intron retention (PSI).
#'
#' NB: SpliceWiz separately considers known "RI" and novel "IR" events 
#'   separately:
#' * For **IR** (all introns are putative IR events), differential analysis
#' is performed on IR-ratio, which is similar to PSI, but the "excluded isoform"
#' includes **all** spliced transcripts that contain an overlapping intron, as
#' estimated via the `SpliceMax` and `SpliceOver` methods - see [collateData].
#' * For **RI** (annotated retained introns), differential analysis is performed
#' on PSI, whereby the "excluded" isoform is one that splice across the 
#' same **exact** intron as the intron of interest. 
#'
#' SpliceWiz considers "included" counts as those that represent abundance of 
#' the "included" isoform, whereas "excluded" counts represent the abundance of 
#' the "excluded" isoform.
#' For consistency, it applies a convention whereby
#' the "included" transcript is one where its splice junctions
#' are by definition shorter than those of "excluded" transcripts.
#' Specifically, this means the included / excluded isoforms are as follows:
#'
#' | EventType | Included | Excluded |
#' | :---: | :---: | :---: |
#' | IR or RI | Intron Retention | Spliced Intron |
#' | MXE | Upstream exon inclusion | Downstream exon inclusion |
#' | SE | Exon inclusion | Exon skipping |
#' | AFE | Downstream exon usage | Upstream exon usage |
#' | ALE | Upstream exon usage | Downstream exon usage |
#' | A5SS | Downstream 5'-SS | Upstream 5'-SS |
#' | A3SS | Upstream 3'-SS | Downstream 3'-SS |
#'
#' @param se The \linkS4class{NxtSE} object created by `makeSE()`. To reduce
#'   runtime and avoid excessive multiple testing, consider filtering
#'   the object using [applyFilters]
#' @param test_factor The condition type which contains the
#'   contrasting variable
#' @param test_nom The nominator condition to test for differential ASE. Usually
#'   the "treatment" condition
#' @param test_denom The denominator condition to test against for differential
#'   ASE. Usually the "control" condition
#' @param batch1,batch2 (Optional, limma and DESeq2 only) One or two condition
#'   types containing batch information to account for.
#' @param filter_antiover,filter_antinear Whether to remove novel IR events that
#'   overlap over or near anti-sense genes. Default will exclude antiover but
#'   not antinear introns. These are ignored if strand-specific RNA-seq 
#'   protocols are used.
#' @param n_threads (DESeq2 only) How many threads to use for DESeq2
#'   based analysis.
#' @return For all methods, a data.table containing the following:
#'   * EventName: The name of the ASE event. This identifies each ASE
#'     in downstream functions including [makeMeanPSI], [makeMatrix],
#'     and [plotCoverage]
#'   * EventType: The type of event. See details section above.
#'   * EventRegion: The genomic coordinates the event occupies. This spans the
#'     most upstream and most downstream splice junction involved in the ASE,
#'     and is use to guide the [plotCoverage] function.
#'   * NMD_direction: Indicates whether one isoform is a NMD substrate. +1 means
#'     included isoform is NMD, -1 means the excluded isoform is NMD, and 0
#'     means there is no change in NMD status (i.e. both / neither are NMD)
#'   * AvgPSI_nom, Avg_PSI_denom: the average percent spliced in / percent
#'     IR levels for the two conditions being contrasted. `nom` and `denom` in
#'     column names are replaced with the condition names
#'
#'   **limma specific output**
#'   * logFC, AveExpr, t, P.Value, adj.P.Val, B: limma topTable columns of
#'     differential ASE. See [limma::topTable] for details.
#'   * inc/exc_(logFC, AveExpr, t, P.Value, adj.P.Val, B): limma results
#'     for differential testing for raw included / excluded counts only
#'
#'   **DESeq2 specific output**
#'   * baseMean, log2FoldChange, lfcSE, stat, pvalue, padj:
#'     DESeq2 results columns for differential ASE; see [DESeq2::results] for
#'     details.
#'   * inc/exc_(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj):
#'     DESeq2 results for differential testing for
#'     raw included / excluded counts only
#'
#'   **satuRn specific output**
#'   * estimates, se, df, t, pval, regular_FDR:
#'     estimated log-odds ratio, standard error, degrees of freedom, (Wald) t
#'     statistic, nominal p-value and associated false discovery rate
#'   * empirical_pval, empirical_FDR: nominal p value and associated FDR
#'     computed by estimating the null distribution of the test statistic
#'     empirically (by satuRn).
#'
#'   **DoubleExp specific output**
#'   * MLE_nom, MLE_denom: Maximum likelihood expectation of PSI values for the 
#"     two groups. `nom` and
#'     `denom` in column names are replaced with the condition names
#'   * MLE_LFC: Log2-fold change of the MLE
#'   * P.Value, adj.P.Val: Nominal and BH-adjusted P values
#'   * n_eff: Number of effective samples (i.e. non-zero or non-unity PSI)
#'   * mDepth: Mean Depth of splice coverage in each of the two groups.
#'   * Dispersion_Reduced, Dispersion_Full: Dispersion values for reduced and
#'     full models. See [DoubleExpSeq::DBGLM1] for details.
#' @examples
#' # Load the NxtSE object and set up the annotations
#' # - see ?makeSE on example code of generating this NxtSE object
#' se <- SpliceWiz_example_NxtSE()
#'
#' colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' colData(se)$replicate <- rep(c("P","Q","R"), 2)
#'
#' require("limma")
#' res_limma <- ASE_limma(se, "treatment", "A", "B")
#'
#' require("DoubleExpSeq")
#' res_DES <- ASE_DoubleExpSeq(se, "treatment", "A", "B")
#'
#' require("satuRn")
#' res_sat <- ASE_satuRn(se, "treatment", "A", "B")
#' 
#' require("DESeq2")
#' res_DESeq <- ASE_DESeq(se, "treatment", "A", "B")
#' 
#' # Time series example
#'
#' colData(se)$timepoint <- rep(c(1,2,3), each = 2)
#' colData(se)$batch <- rep(c("1", "2"), 3)
#' res_DESeq_timeseries <- ASE_DESeq(se, "timepoint")
#' 
#' @name ASE-methods
#' @references
#' Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015).
#' 'limma powers differential expression analyses for RNA-sequencing and
#' microarray studies.' Nucleic Acids Research, 43(7), e47.
#' \url{https://doi.org/10.1093/nar/gkv007}
#'
#' Love MI, Huber W, Anders S (2014). 'Moderated estimation of fold change and
#' dispersion for RNA-seq data with DESeq2.' Genome Biology, 15, 550.
#' \url{https://doi.org/10.1186/s13059-014-0550-8}
#'
#' Ruddy S, Johnson M, Purdom E (2016). 'Shrinkage of dispersion parameters in
#' the binomial family, with application to differential exon skipping.'
#' Ann. Appl. Stat. 10(2): 690-725.
#' \url{https://doi.org/10.1214/15-AOAS871}
#'
#' Gilis J, Vitting-Seerup K, Van den Berge K, Clement L (2021). 'Scalable
#' analysis of differential transcript usage for bulk and single-cell
#' RNA-sequencing applications.' F1000Research 2021, 10:374.
#' \url{https://doi.org/10.12688/f1000research.51749.1}
#' @md
NULL

#' @describeIn ASE-methods Use limma to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
ASE_limma <- function(se, test_factor, test_nom, test_denom,
        batch1 = "", batch2 = "",
        filter_antiover = TRUE, filter_antinear = FALSE) {

    .check_package_installed("limma", "3.44.0")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, batch1, batch2)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear)

    .log("Performing limma contrast for included / excluded counts separately",
        "message")
    res.limma2 <- .ASE_limma_contrast(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2)
    res.inc <- res.limma2[grepl(".Included", get("EventName"))]
    res.inc[, c("EventName") :=
        sub(".Included","",get("EventName"), fixed=TRUE)]
    res.inc <- res.inc[get("AveExpr") > 1]   # Filter as 0/5 is not diff to 0/10
    res.exc <- res.limma2[grepl(".Excluded", get("EventName"))]
    res.exc[, c("EventName") :=
        sub(".Excluded","",get("EventName"), fixed=TRUE)]
    res.exc <- res.exc[get("AveExpr") > 1]

    .log("Performing limma contrast for included / excluded counts together",
        "message")
    rowData <- as.data.frame(rowData(se_use))
    se_use <- se_use[rowData$EventName %in% res.inc$EventName &
        rowData$EventName %in% res.exc$EventName,]
    res.ASE <- .ASE_limma_contrast_ASE(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2)
    res.ASE[res.inc, on = "EventName",
        paste("Inc", colnames(res.inc)[seq_len(6)], sep=".") :=
        list(get("i.logFC"), get("i.AveExpr"), get("i.t"),
            get("i.P.Value"), get("i.adj.P.Val"), get("i.B"))]
    res.ASE[res.exc, on = "EventName",
        paste("Exc", colnames(res.exc)[seq_len(6)], sep=".") :=
        list(get("i.logFC"), get("i.AveExpr"), get("i.t"),
            get("i.P.Value"), get("i.adj.P.Val"), get("i.B"))]
    setorderv(res.ASE, "B", order = -1)
    res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, 
        test_nom, test_denom)
    return(res.ASE)
}

#' @describeIn ASE-methods Use DESeq2 to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
ASE_DESeq <- function(se, test_factor, test_nom, test_denom,
        batch1 = "", batch2 = "",
        n_threads = 1,
        filter_antiover = TRUE, filter_antinear = FALSE) {
    .check_package_installed("DESeq2", "1.30.0")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, batch1, batch2,
        allowTimeSeries = TRUE)
    BPPARAM_mod <- .validate_threads(n_threads)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear)

    .log("Performing DESeq2 contrast for included / excluded counts separately",
        "message")
    res.IncExc <- .ASE_DESeq2_contrast(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM_mod)
    res.inc <- res.IncExc[grepl(".Included", get("EventName"))]
    res.inc[, c("EventName") :=
        sub(".Included","",get("EventName"), fixed=TRUE)]
    res.exc <- res.IncExc[grepl(".Excluded", get("EventName"))]
    res.exc[, c("EventName") :=
        sub(".Excluded","",get("EventName"), fixed=TRUE)]

    .log("Performing DESeq2 contrast for included / excluded counts separately",
        "message")
    rowData <- as.data.frame(rowData(se_use))
    se_use <- se_use[rowData$EventName %in% res.inc$EventName &
        rowData$EventName %in% res.exc$EventName,]
    res.ASE <- .ASE_DESeq2_contrast_ASE(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM_mod)
    res.ASE[res.inc, on = "EventName",
        paste("Inc", colnames(res.inc)[seq_len(6)], sep=".") :=
        list(get("i.baseMean"), get("i.log2FoldChange"), get("i.lfcSE"),
            get("i.stat"), get("i.pvalue"), get("i.padj"))]
    res.ASE[res.exc, on = "EventName",
        paste("Exc", colnames(res.exc)[seq_len(6)], sep=".") :=
        list(get("i.baseMean"), get("i.log2FoldChange"), get("i.lfcSE"),
            get("i.stat"), get("i.pvalue"), get("i.padj"))]
    res.ASE <- res.ASE[!is.na(get("pvalue"))]
    setorder(res.ASE, "pvalue")
    if(is_valid(test_nom)) {
        res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, 
            test_nom, test_denom)
    } else {
        condlist <- as.list(sort(unique(
            unlist(colData(se)[, test_factor]
        ))))
        if(length(condlist) > 6) condlist <- condlist[seq_len(6)]
        res.ASE <- .ASE_add_diag_multi(res.ASE, se_use, test_factor, condlist)
    }
    return(res.ASE)
}

#' @describeIn ASE-methods Use DoubleExpSeq to perform differential ASE analysis
#'   of a filtered NxtSE object (uses double exponential beta-binomial model)
#'   to estimate group dispersions, followed by LRT
#' @export
ASE_DoubleExpSeq <- function(se, test_factor, test_nom, test_denom,
        # batch1 = "", batch2 = "",
        filter_antiover = TRUE, filter_antinear = FALSE) {

    .check_package_installed("DoubleExpSeq", "1.1")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, "", "")
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear)

    .log("Running DoubleExpSeq::DBGLM1() on given data", "message")
    res.ASE <- .ASE_DoubleExpSeq_contrast_ASE(se_use,
        test_factor, test_nom, test_denom)

    res.cols <- c(
        paste("MLE", test_nom, sep="_"), paste("MLE", test_denom, sep="_"),
        "P.Value", "adj.P.Val", "n_eff",
        paste("mDepth", test_nom, sep="_"),
        paste("mDepth", test_denom, sep="_"),
        "Dispersion_Reduced", "Dispersion_Full"
    )
    colnames(res.ASE)[-1] <- res.cols

    res.ASE[, c("MLE_LFC") := (
        qlogis(res.ASE[,get(paste("MLE", test_nom, sep="_"))]) -
        qlogis(res.ASE[,get(paste("MLE", test_denom, sep="_"))])
    ) / log(2)]

    res.ASE <- res.ASE[, c("EventName", res.cols[c(1,2)], "MLE_LFC",
        res.cols[seq(3,9)]), with = FALSE]

    res.ASE <- res.ASE[!is.na(get("P.Value"))]
    setorderv(res.ASE, "P.Value")
    res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, 
        test_nom, test_denom)
    return(res.ASE)
}

#' @describeIn ASE-methods Use satuRn to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
ASE_satuRn <- function(se, test_factor, test_nom, test_denom,
        batch1 = "", batch2 = "",
        filter_antiover = TRUE, filter_antinear = FALSE) {

    .check_package_installed("satuRn", "1.4.2")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, batch1, batch2)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear)

    .log("Performing satuRn contrast for included / excluded counts",
        "message")
    rowData <- as.data.frame(rowData(se_use))
    res.ASE <- .ASE_satuRn_contrast(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2)
    res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, 
        test_nom, test_denom)
    return(res.ASE)
}

############################## INTERNALS #######################################
# helper functions:

# Check arguments are valid
.ASE_check_args <- function(colData, test_factor,
        test_nom, test_denom, batch1, batch2, 
        allowTimeSeries = FALSE
) {
    if(!is_valid(test_factor)) {
        .log("test_factor must be defined")
    } else if (!allowTimeSeries & 
            (!is_valid(test_nom) | !is_valid(test_denom))) {
        .log("test_nom, test_denom must be defined")
    } else if(allowTimeSeries & !is_valid(test_nom)) {
        if(!is.numeric(colData[, test_factor])) {
            .log(paste(
                test_factor, "must be numeric for time series analysis"
            ))
        }
    }
    if(!(test_factor %in% colnames(colData))) {
        .log("test_factor is not a condition in colData")
    }
    if(is_valid(test_nom) && !any(colData[, test_factor] == test_nom)) {
        .log("test_nom is not found in any samples")
    }
    if(is_valid(test_denom) && !any(colData[, test_factor] == test_denom)) {
        .log("test_denom is not found in any samples")
    }
    if(batch1 != "") {
        if(!(batch1 %in% colnames(colData))) {
            .log("batch1 is not a condition in colData")
        }
        if(test_factor == batch1) {
            .log("batch1 and test_factor are the same")
        }
    }
    if(batch2 != "") {
        if(!(batch2 %in% colnames(colData))) {
            .log("batch2 is not a condition in colData")
        }
        if(test_factor == batch2) {
            .log("batch2 and test_factor are the same")
        }
    }
    if(batch1 != "" & batch2 != "") {
        if(batch1 == batch2) {
            .log("batch1 and batch2 are the same")
        }
    }
    return(TRUE)
}

# Filter antiover and antinear
.ASE_filter <- function(se, filter_antiover, filter_antinear) {
    se_use <- se
    if(filter_antiover) {
        se_use <- se_use[
            !grepl("anti-over", rowData(se_use)$EventName),]
    }
    if(filter_antinear) {
        se_use <- se_use[
            !grepl("anti-near", rowData(se_use)$EventName),]
    }
    se_use <- se_use[
        !grepl("known-exon", rowData(se_use)$EventName),]
    return(se_use)
}

.ASE_limma_contrast <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
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

    condition_factor <- factor(colData[, test_factor])
    if(batch2 != "") {
        batch2_factor <- colData[, batch2]
        batch1_factor <- colData[, batch1]
        design1 <- model.matrix(~0 + batch1_factor + batch2_factor +
            condition_factor)
    } else if(batch1 != "") {
        batch1_factor <- colData[, batch1]
        design1 <- model.matrix(~0 + batch1_factor + condition_factor)
    } else {
        design1 <- model.matrix(~0 + condition_factor)
    }
    contrast <- rep(0, ncol(design1))
    contrast_a <- paste0("condition_factor", test_nom)
    contrast_b <- paste0("condition_factor", test_denom)
    contrast[which(colnames(design1) == contrast_b)] <- -1
    contrast[which(colnames(design1) == contrast_a)] <- 1

    countData_use <- limma::voom(countData, design1)
    fit <- limma::lmFit(countData_use$E, design = design1)
    fit <- limma::contrasts.fit(fit, contrast)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, number = nrow(countData_use$E))
    res$EventName <- rownames(res)

    res$AveExpr <- res$AveExpr - min(res$AveExpr)
    res <- as.data.table(res)
    
    rm(fit, countData, countData_use)
    gc()
    return(res)
}

.ASE_limma_contrast_ASE <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
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

    condition_factor <- factor(colData[, test_factor])
    ASE <- colData[, "ASE"]
    if(batch2 != "") {
        batch2_factor <- colData[, batch2]
        batch1_factor <- colData[, batch1]
        design1 <- model.matrix(~0 + batch1_factor + batch2_factor +
            condition_factor + condition_factor:ASE)
    } else if(batch1 != "") {
        batch1_factor <- colData[, batch1]
        design1 <- model.matrix(~0 + batch1_factor + condition_factor +
            condition_factor:ASE)
    } else {
        design1 <- model.matrix(~0 + condition_factor + condition_factor:ASE)
    }
    colnames(design1) <- sub(":",".",colnames(design1))
    contrast <- rep(0, ncol(design1))
    contrast_a <- paste0("condition_factor", test_nom, ".ASEIncluded")
    contrast_b <- paste0("condition_factor", test_denom, ".ASEIncluded")
    contrast[which(colnames(design1) == contrast_b)] <- -1
    contrast[which(colnames(design1) == contrast_a)] <- 1

    countData_use <- limma::voom(countData, design1, lib.size = 1)

    fit <- limma::lmFit(countData_use$E, design = design1)

    fit <- limma::contrasts.fit(fit, contrast)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, number = nrow(countData_use$E))
    res$EventName <- rownames(res)
    res$AveExpr <- res$AveExpr - min(res$AveExpr)
    res <- as.data.table(res)
    
    rm(fit, countData, countData_use)
    gc()
    return(res)
}

.ASE_DESeq2_contrast <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM) {
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
    if(batch2 != "") {
        dds_formula <- paste0("~", paste(
            batch1, batch2, test_factor,
            sep="+"))
        dds_formula_reduced <- paste0("~", paste(
            batch1, batch2, sep="+"))
    } else if(batch1 != "") {
        dds_formula <- paste0("~", paste(
            batch1, test_factor,
            sep="+"))
        dds_formula_reduced <- paste0("~", paste(
            batch1, sep="+"))
    } else {
        dds_formula <- paste0("~", test_factor)
        dds_formula_reduced <- paste0("~1")
    }

    mode(countData) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(countData),
        colData = .DESeq_colData_sanitise(colData),
        design = as.formula(dds_formula)
    )
    message("ASE_DESeq: Profiling expression of Included and Excluded counts")

    if(!is_valid(test_nom) | !is_valid(test_denom)) {
        dds <- DESeq2::DESeq(dds, test = "LRT", 
            reduced = as.formula(dds_formula_reduced),
            parallel = TRUE, BPPARAM = BPPARAM)
        res <- as.data.frame(DESeq2::results(dds,
            parallel = TRUE, BPPARAM = BPPARAM)
        )        
    } else {
        dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)
        res <- as.data.frame(DESeq2::results(dds,
            contrast = c(test_factor, test_nom, test_denom),
            parallel = TRUE, BPPARAM = BPPARAM)
        )    
    }
    res$EventName <- rownames(res)

    rm(dds, countData)
    gc()
    return(as.data.table(res))
}

.ASE_DESeq2_contrast_ASE <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM) {
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

    if(batch2 != "") {
        dds_formula <- paste0("~", paste(
            batch1, batch2, test_factor,
                paste0(test_factor, ":ASE"),
            sep="+"))
        dds_formula_reduced <- paste0("~", paste(
            batch1, batch2, test_factor,
                # paste0(test_factor, ":ASE"),
            sep="+"))
    } else if(batch1 != "") {
        dds_formula <- paste0("~", paste(
            batch1, test_factor,
            paste0(test_factor, ":ASE"),
            sep="+"))
        dds_formula_reduced <- paste0("~", paste(
            batch1, test_factor,
            # paste0(test_factor, ":ASE"),
            sep="+"))
    } else {
        dds_formula <- paste0("~", paste(
            test_factor,
            paste0(test_factor, ":ASE"),
            sep="+"))
        dds_formula_reduced <- paste0("~", paste(
            test_factor,
            # paste0(test_factor, ":ASE"),
            sep="+"))
    }

    mode(countData) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = countData,
        colData = .DESeq_colData_sanitise(colData),
        design = as.formula(dds_formula)
    )
    DESeq2::sizeFactors(dds) <- 1
    message("ASE_DESeq: Profiling differential ASE")
    if(!is_valid(test_nom) | !is_valid(test_denom)) {
        dds <- DESeq2::DESeq(dds, test = "LRT",
            reduced = as.formula(dds_formula_reduced),
            parallel = TRUE, BPPARAM = BPPARAM)
        res <- as.data.frame(DESeq2::results(dds,
            parallel = TRUE, BPPARAM = BPPARAM)
        )
    } else {
        dds <- DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)
        res <- as.data.frame(DESeq2::results(dds,
            list(
                paste0(test_factor, test_nom, ".ASEIncluded"),
                paste0(test_factor, test_denom, ".ASEIncluded")
            ),
            parallel = TRUE, BPPARAM = BPPARAM)
        )
    }
    res$EventName <- rownames(res)
    
    rm(dds, countData)
    gc()
    return(as.data.table(res))
}

.ASE_DoubleExpSeq_contrast_ASE <- function(se, test_factor,
    test_nom, test_denom) {

    # NB add pseudocounts
    pseudocount <- 1
    y <- as.matrix(assay(se, "Included") + pseudocount)
    m <- as.matrix(
        assay(se, "Included") + assay(se, "Excluded") + 2 * pseudocount
    )
    colData <- as.data.frame(colData(se))
    groups <- factor(colData[, test_factor])
    shrink.method <- "WEB"

    # carry over error from DoubleExpSeq version 1.1
    contrast.first <- which(unique(groups) == test_nom)
    contrast.second <- which(unique(groups) == test_denom)

    res <- DoubleExpSeq::DBGLM1(
        as.matrix(y), as.matrix(m), groups, shrink.method,
        contrast=c(contrast.first,contrast.second),
        fdr.level=0.05, use.all.groups=TRUE)

    rm(y, m)
    gc()
    return(cbind(data.table(EventName = rownames(res$All)),
        as.data.table(res$All)))
}

.ASE_satuRn_contrast <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
    countData <- as.matrix(rbind(assay(se, "Included"),
        assay(se, "Excluded")))
    rowData <- as.data.frame(rowData(se))
    colData <- colData(se)
    rownames(colData) <- colnames(se)
    colnames(countData) <- rownames(colData)
    
    txInfo <- as.data.frame(matrix(data = NA, nrow = 2 * nrow(se), ncol = 2))
    colnames(txInfo) <- c("isoform_id", "gene_id")
    txInfo$isoform_id <- c(
        paste("Inc", rowData$EventName, sep=":"), 
        paste("Exc", rowData$EventName, sep=":")
    )
    txInfo$gene_id <- rep(rowData$EventName, 2)
    
    rownames(countData) <- txInfo$isoform_id

    condition_factor <- factor(colData[, test_factor])
    if(batch2 != "") {
        batch2_factor <- colData[, batch2]
        batch1_factor <- colData[, batch1]
        formula1 <- ~0 + batch1_factor + batch2_factor +
            condition_factor
    } else if(batch1 != "") {
        batch1_factor <- colData[, batch1]
        formula1 <- ~0 + batch1_factor + condition_factor
    } else {
        formula1 <- ~0 + condition_factor
    }
    design1 <- model.matrix(formula1)
    contrast <- rep(0, ncol(design1))
    contrast_a <- paste0("condition_factor", test_nom)
    contrast_b <- paste0("condition_factor", test_denom)
    contrast[which(colnames(design1) == contrast_b)] <- -1
    contrast[which(colnames(design1) == contrast_a)] <- 1

    contrast <- matrix(contrast, ncol = 1)
    colnames(contrast) <- "Contrast1"

    sumExp <- SummarizedExperiment(
        assays = list(counts = countData),
        colData = colData,
        rowData = txInfo
    )
    sumExp <- satuRn::fitDTU(
        object = sumExp,
        formula = formula1,
        verbose = FALSE
    )
    sumExp <- satuRn::testDTU(
        object = sumExp,
        contrasts = contrast,
        diagplot1 = FALSE,
        diagplot2 = FALSE,
        sort = FALSE
    )

    res <- as.data.frame(rowData(sumExp)[["fitDTUResult_Contrast1"]])
    res <- res[substr(rownames(res), 1, 4) == "Inc:",]
    rownames(res) <- substr(rownames(res), 5, nchar(rownames(res)))
    res.ASE <- as.data.table(cbind(
        data.frame(EventName = rownames(res), stringsAsFactors = FALSE),
        res
    ))
    setorderv(res.ASE, "pval")
    res.ASE <- res.ASE[!is.na(get("estimates"))]
    res.ASE[, c("log2estimates") := list(get("estimates") / log(2))]

    rm(res, sumExp, countData)
    gc()
    return(res.ASE)
}

.ASE_add_diag <- function(res, se, test_factor, test_nom, test_denom) {
    rowData <- as.data.frame(rowData(se))
    rowData.DT <- as.data.table(rowData[,
        c("EventName","EventType","EventRegion", "NMD_direction")])
    diag <- makeMeanPSI(se, res$EventName,
        test_factor, list(test_nom, test_denom))
    colnames(diag)[2:3] <- c(paste0("AvgPSI_", test_nom),
        paste0("AvgPSI_", test_denom))
    res <- cbind(
        res[,c("EventName")],
        as.data.table(round(diag[,-1], 4)),
        res[,-c("EventName")]
    )
    res <- rowData.DT[res, on = "EventName"]
    return(res)
}

.ASE_add_diag_multi <- function(
        res, se, test_factor, 
        conditionList
) {
    rowData <- as.data.frame(rowData(se))
    rowData.DT <- as.data.table(rowData[,
        c("EventName","EventType","EventRegion", "NMD_direction")])
    diag <- makeMeanPSI(se, res$EventName,
        test_factor, conditionList)
    for(i in seq_len(length(conditionList))) {
        colnames(diag)[i+1] <- paste0("AvgPSI_", conditionList[i])
    }
    res <- cbind(
        res[,c("EventName")],
        as.data.table(round(diag[,-1], 4)),
        res[,-c("EventName")]
    )
    res <- rowData.DT[res, on = "EventName"]
    return(res)
}

.DESeq_colData_sanitise <- function(colData) {
    for(i in seq_len(ncol(colData))) {
        if(is(colData[,i], "character")) {
            colData[, i] <- factor(colData[, i])
        }
    }
    colData
}
