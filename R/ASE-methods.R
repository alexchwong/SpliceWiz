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
#' Time series are supported by SpliceWiz to a limited extent. 
#' Time series analysis can be performed via limma or DESeq2.
#' For limma time-series analysis, use `ASE_limma_timeseries()`, specifying
#' the `test_factor` as the column of numeric values containing 
#' time series data. For DESeq, time series differential analysis can be 
#' activated using the `ASE_DESeq()` function, again specifying `test_factor`
#' as the column containing time series data (and leaving `test_nom`
#' and `test_denom` parameters blank). See examples below.
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
#' NB: SpliceWiz measures intron retention events using two different 
#'   approaches, the choice of which is left to the user - see [ASE-methods]:
#' * **IR** (intron retention) events: considers all introns to be potentially
#' retained. Given in most scenarios there may be uncertainty as to which of the
#' many mutually-overlapping introns are spliced to produce the major isoform,
#' SpliceWiz adopts the IRFinder approach by using the IR-ratio. The "included"
#' isoform is the relative abundance of the IR-transcript, as approximated by
#' the trimmed-mean depth of coverage across the intron (excluding outliers
#' including exons of other transcripts, intronic elements such as snoRNAs, 
#' etc). The "excluded isoform" includes **all** spliced transcripts that 
#' contain an overlapping intron, as estimated via SpliceWiz's `SpliceOver` and 
#' IRFinder's `SpliceMax` methods - see [collateData].
#' * **RI** (annotated retained introns) considers only annotated retained 
#' introns, i.e., those annotated within the given reference. These are
#' quantified using PSI, considering the included (IR-transcript) and excluded
#' (splicing of the exact intron) as binary alternatives. 
#'
#' SpliceWiz considers "included" counts as those that represent abundance of 
#' the "included" isoform, whereas "excluded" counts represent the abundance of 
#' the "excluded" isoform.
#' To allow comparison between modalities, SpliceWiz applies a convention 
#' whereby the "included" transcript is one where its splice junctions
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
#' @param test_factor The column name in the sample annotation `colData(se)` 
#'   that contains the desired variables to be contrasted. For 
#'   `ASE_limma_timeseries()` and `ASE_DESeq()` (when `test_nom` and 
#'   `test_denom` parameters are left blank), `test_factor` must contain 
#'   numerical values representing the time variable.
#' @param test_nom The nominator condition to test for differential ASE. Usually
#'   the "treatment" condition
#' @param test_denom The denominator condition to test against for differential
#'   ASE. Usually the "control" condition
#' @param batch1,batch2 (Optional, limma and DESeq2 only) One or two condition
#'   types containing batch information to account for.
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
#' @param n_threads (DESeq2 only) How many threads to use for DESeq2
#'   based analysis.
#' @param degrees_of_freedom (default `1`) The complexity of time series trends
#'   modeled by `ASE_limma_timeseries`. E.g., `1` will only model linear trends,
#'   whereas `2` extends the capacity for quadratic trends, `3` for cubic
#'   trends, etc.
#' @param filterByMinCPM (default `0`) In `ASE_satuRn()`, Included/Excluded
#'   counts will be filtered using this value as the threshold prior to satuRn
#'   analysis. Filtering is performed using `edgeR::filterByExpr()` parsing
#'   this parameter into its `min.count` parameter.
#' @return For all methods, a data.table containing the following:
#'   * EventName: The name of the ASE event. This identifies each ASE
#'     in downstream functions including [makeMeanPSI], [makeMatrix],
#'     and [plotCoverage]
#'   * EventType: The type of event. See details section above.
#'   * EventRegion: The genomic coordinates the event occupies. This spans the
#'     most upstream and most downstream splice junction involved in the ASE,
#'     and is use to guide the [plotCoverage] function.
#'   * flags: Indicates which isoforms are NMD substrates and/or which are
#'     formed by novel splicing only.
#'   * AvgPSI_nom, Avg_PSI_denom: the average percent spliced in / percent
#'     IR levels for the two conditions being contrasted. `nom` and `denom` in
#'     column names are replaced with the condition names. Note this is a
#'     geometric mean, based on the arithmetic mean of logit PSI values.
#'   * deltaPSI: The difference in PSI between the mean values of the two
#'     conditions.
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
#' # Limma analysis (counts modeled using log-normal distribution)
#'
#' require("limma")
#' res_limma <- ASE_limma(se, "treatment", "A", "B")
#'
#' # DoubleExpSeq analysis (counts modeled using beta binomial distribution)
#'
#' require("DoubleExpSeq")
#' res_DES <- ASE_DoubleExpSeq(se, "treatment", "A", "B")
#'
#' # satuRn analysis (quasi binomial), 
#' #   Filtering counts using 1 count per million as threshold prior to analysis
#'
#' require("satuRn")
#' require("edgeR")
#'
#' res_sat <- ASE_satuRn(se, "treatment", "A", "B", filterByMinCPM = 1)
#' 
#' # DESeq2 analysis (counts modeled using negative binomial distribution)
#'
#' require("DESeq2")
#' res_DESeq <- ASE_DESeq(se, "treatment", "A", "B")
#' 
#' # Time series examples
#'
#' colData(se)$timepoint <- rep(c(1,2,3), each = 2)
#' colData(se)$batch <- rep(c("1", "2"), 3)
#'
#' res_limma_timeseries <- ASE_limma_timeseries(se, "timepoint")
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
#' See also [ASE-methods] for explanations of results output.
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
#' @return For all methods, a data.table containing the following:
#'   * EventName: The name of the ASE event. This identifies each ASE
#'     in downstream functions including [makeMeanPSI], [makeMatrix],
#'     and [plotCoverage]
#'   * EventType: The type of event. See details section above.
#'   * EventRegion: The genomic coordinates the event occupies. This spans the
#'     most upstream and most downstream splice junction involved in the ASE,
#'     and is use to guide the [plotCoverage] function.
#'   * flags: Indicates which isoforms are NMD substrates and/or which are
#'     formed by novel splicing only.
#'   * AvgPSI_nom, Avg_PSI_denom: the average percent spliced in / percent
#'     IR levels for the two conditions being contrasted. `nom` and `denom` in
#'     column names are replaced with the condition names. Note this is a
#'     geometric mean, based on the arithmetic mean of logit PSI values.
#'   * deltaPSI: The difference in PSI between the mean values of the two
#'     conditions.
#'
#'   **edgeR specific output**
#'   * logFC, AveExpr, t, P.Value, adj.P.Val, B: limma topTable columns of
#'     differential ASE. See [limma::topTable] for details.
#'   * inc/exc_(logFC, AveExpr, t, P.Value, adj.P.Val, B): limma results
#'     for differential testing for raw included / excluded counts only
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

#' @describeIn ASE-methods Use limma to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
ASE_limma <- function(se, test_factor, test_nom, test_denom,
        batch1 = "", batch2 = "",
        IRmode = c("all", "annotated", "annotated_binary"),
        filter_antiover = TRUE, filter_antinear = FALSE) {

    .check_package_installed("limma", "3.44.0")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, batch1, batch2)
    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)

    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")

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

#' @describeIn ASE-methods Use edgeR to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
ASE_edgeR <- function(se, test_factor, test_nom, test_denom,
        batch1 = "", batch2 = "",
        IRmode = c("all", "annotated", "annotated_binary"),
        filter_antiover = TRUE, filter_antinear = FALSE) {

    .check_package_installed("edgeR", "3.32.0")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, batch1, batch2)
    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)

    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")

    .log("Performing edgeR contrast for included / excluded counts separately",
        "message")
    res.edgeR <- .ASE_edgeR_contrast(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2)
    res.inc <- res.edgeR[grepl(".Included", get("EventName"))]
    res.inc[, c("EventName") :=
        sub(".Included","",get("EventName"), fixed=TRUE)]
    # res.inc <- res.inc[get("AveExpr") > 1]   # Filter as 0/5 is not diff to 0/10
    res.exc <- res.edgeR[grepl(".Excluded", get("EventName"))]
    res.exc[, c("EventName") :=
        sub(".Excluded","",get("EventName"), fixed=TRUE)]
    # res.exc <- res.exc[get("AveExpr") > 1]

    .log("Performing edgeR contrast for included / excluded counts together",
        "message")
    rowData <- as.data.frame(rowData(se_use))
    se_use <- se_use[rowData$EventName %in% res.inc$EventName &
        rowData$EventName %in% res.exc$EventName,]
    res.ASE <- .ASE_edgeR_contrast_ASE(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2)
    res.ASE[res.inc, on = "EventName",
        paste("Inc", colnames(res.inc)[seq_len(5)], sep=".") :=
        list(get("i.logFC"), get("i.logCPM"), get("i.F"),
            get("i.PValue"), get("i.FDR"))]
    res.ASE[res.exc, on = "EventName",
        paste("Exc", colnames(res.inc)[seq_len(5)], sep=".") :=
        list(get("i.logFC"), get("i.logCPM"), get("i.F"),
            get("i.PValue"), get("i.FDR"))]
    setorderv(res.ASE, "F", order = -1)
    res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, 
        test_nom, test_denom)
    return(res.ASE)
}

#' @describeIn ASE-GLM-edgeR Use edgeR to fit counts and ASE models with a
#'   given design formula
#' @export
fitASE_edgeR <- function(
        se,
        strModelFormula, strASEFormula,
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
    fit_ASE <- .ASE_edgeR_fitASE(se_use, models$design_ASE)

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
    fit_ASE <- .ASE_edgeR_fitASE(se_use, model_ASE)

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
    coef_IncExc = ncol(model_IncExc), contrast_IncExc = NULL,
    coef_ASE = ncol(model_ASE), contrast_ASE = NULL
) {
    se_use <- se[rownames(se) %in% rownames(fit$ASE)]
    qlf_IncExc <- edgeR::glmQLFTest(
        fit$IncExc,
        coef = coef_IncExc, contrast = contrast_IncExc
    )
    res_IncExc <- edgeR::topTags(qlf_IncExc, n = nrow(fit$fit_IncExc))
    res_IncExc$table$EventName <- rownames(res_IncExc$table)

    res_IncExc <- as.data.table(res_IncExc$table)
    res.inc <- res_IncExc[grepl(".Included", get("EventName"))]
    res.inc[, c("EventName") :=
        sub(".Included","",get("EventName"), fixed=TRUE)]
    res.exc <- res_IncExc[grepl(".Excluded", get("EventName"))]
    res.exc[, c("EventName") :=
        sub(".Excluded","",get("EventName"), fixed=TRUE)]
        
    qlf_ASE <- edgeR::glmQLFTest(
        fit$ASE,
        coef = coef_ASE, contrast = contrast_ASE
    )
    res_ASE <- edgeR::topTags(qlf_ASE, n = nrow(fit$fit_IncExc))
    res_ASE$table$EventName <- rownames(res_ASE)
    
    res_ASE <- as.data.table(res_ASE$table)
    colnames(res.inc)[-ncol(res.inc)] <- paste(
        "Inc", colnames(res.inc)[-ncol(res.inc)], sep = ".")
    colnames(res.exc)[-ncol(res.exc)] <- paste(
        "Inc", colnames(res.exc)[-ncol(res.exc)], sep = ".")
    
    res_ASE <- res_ASE[res.inc, on = "EventName"]
    res_ASE <- res_ASE[res.exc, on = "EventName"]
    setorderv(res_ASE, "F", order = -1)
    
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

#' @describeIn ASE-methods Use limma to perform differential ASE analysis of
#'   a filtered NxtSE object (time series)
#' @export
ASE_limma_timeseries <- function(se, test_factor,
        batch1 = "", batch2 = "",
        degrees_of_freedom = 1,
        IRmode = c("all", "annotated", "annotated_binary"),
        filter_antiover = TRUE, filter_antinear = FALSE
) {

    .check_package_installed("limma", "3.44.0")
    .ASE_check_args(colData(se), test_factor, "", "",
        batch1, batch2, allowTimeSeries = TRUE)
    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)

    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")

    .log("Performing limma contrast for included / excluded counts separately",
        "message")
    res.limma2 <- .ASE_limma_contrast_ts(se_use,
        test_factor,
        batch1, batch2, degrees_of_freedom)
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
    res.ASE <- .ASE_limma_contrast_ASE_ts(se_use,
        test_factor, # test_nom, test_denom,
        batch1, batch2, degrees_of_freedom)
    cols <- seq(
        ncol(res.inc) - 4,
        ncol(res.inc) - 1        
    )
    if(degrees_of_freedom == 1) {
        res.ASE[res.inc, on = "EventName",
            paste("Inc", colnames(res.inc)[seq_len(6)], sep=".") :=
            list(get("i.logFC"), get("i.AveExpr"), get("i.t"),
                get("i.P.Value"), get("i.adj.P.Val"), get("i.B"))]
        res.ASE[res.exc, on = "EventName",
            paste("Exc", colnames(res.exc)[seq_len(6)], sep=".") :=
            list(get("i.logFC"), get("i.AveExpr"), get("i.t"),
                get("i.P.Value"), get("i.adj.P.Val"), get("i.B"))]
    
    } else {
        res.ASE[res.inc, on = "EventName",
            paste("Inc", colnames(res.inc)[cols], sep=".") :=
            list(get("i.AveExpr"), get("i.F"),
                get("i.P.Value"), get("i.adj.P.Val"))]
        res.ASE[res.exc, on = "EventName",
            paste("Exc", colnames(res.exc)[cols], sep=".") :=
            list(get("i.AveExpr"), get("i.F"),
                get("i.P.Value"), get("i.adj.P.Val"))]    
    }
    setorderv(res.ASE, "B", order = -1)

    condlist <- as.list(sort(unique(
        unlist(colData(se)[, test_factor]
    ))))
    if(length(condlist) > 6) condlist <- condlist[seq_len(6)]
    res.ASE <- .ASE_add_diag_multi(res.ASE, se_use, test_factor, condlist)

    return(res.ASE)
}


#' @describeIn ASE-methods Use DESeq2 to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
ASE_DESeq <- function(se, test_factor, test_nom, test_denom,
        batch1 = "", batch2 = "",
        n_threads = 1,
        IRmode = c("all", "annotated", "annotated_binary"),
        filter_antiover = TRUE, filter_antinear = FALSE) {
    .check_package_installed("DESeq2", "1.30.0")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, batch1, batch2,
        allowTimeSeries = TRUE)
    BPPARAM_mod <- .validate_threads(n_threads)
    
    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)
    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")

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
        IRmode = c("all", "annotated", "annotated_binary"),
        filter_antiover = TRUE, filter_antinear = FALSE) {

    .check_package_installed("DoubleExpSeq", "1.1")
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, "", "")
        
    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)
    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")

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
        n_threads = 1,
        IRmode = c("all", "annotated", "annotated_binary"),
        filter_antiover = TRUE, filter_antinear = FALSE,
        filterByMinCPM = 0
        # ...
) {

    .check_package_installed("satuRn", "1.4.2")
    # .check_package_installed("edgeR", "3.28.1")
    
    .ASE_check_args(colData(se), test_factor,
        test_nom, test_denom, batch1, batch2)
    BPPARAM_mod <- .validate_threads(n_threads)

    IRmode <- match.arg(IRmode)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, IRmode)

    # Further filtering step using filterByExpr
    if(filterByMinCPM > 0) {
        .check_package_installed("edgeR", "3.28.1")
        countData <- as.matrix(cbind(assay(se_use, "Included"),
            assay(se_use, "Excluded")))
        se_use <- se_use[edgeR::filterByExpr(
            countData, min.count = filterByMinCPM
        ),]
    } else {
        countData <- as.matrix(cbind(assay(se_use, "Included"),
            assay(se_use, "Excluded")))
        if(any(rowSums(countData) == 0)) {
            .log(paste(
                "Events (rows) with all zero counts are removed from analysis."
            ), "warning")
            se_use <- se_use[rowSums(countData) > 0,]
            countData <- as.matrix(cbind(assay(se_use, "Included"),
                assay(se_use, "Excluded")))
        }
    }
    
    if(nrow(se_use) == 0)
        .log("No events for ASE analysis after filtering")
    
    .log("Performing satuRn contrast for included / excluded counts",
        "message")
    rowData <- as.data.frame(rowData(se_use))
    res.ASE <- .ASE_satuRn_contrast(se_use,
        test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM_mod)
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

# Check arguments are valid
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

# Filter antiover and antinear
.ASE_filter <- function(se, filter_antiover, filter_antinear, IRmode) {
    se_use <- se
    if(filter_antiover) {
        se_use <- se_use[
            !grepl("anti-over", rowData(se_use)$EventName),]
    }
    if(filter_antinear) {
        se_use <- se_use[
            !grepl("anti-near", rowData(se_use)$EventName),]
    }
    
    if(IRmode == "all") {
        se_use <- se_use[rowData(se_use)$EventType != "RI",]
    } else if(IRmode == "annotated") {
        se_use <- se_use[rowData(se_use)$EventType != "RI",]
        se_use <- se_use[
            rowData(se_use)$EventType != "IR" | 
            rowData(se_use)$is_annotated_IR,]
    } else {
        se_use <- se_use[rowData(se_use)$EventType != "IR",]    
    }

    return(se_use)
}

.ASE_limma_contrast <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
    in_data <- .ASE_contrast_expr(se, test_factor, 
        test_nom, test_denom,
        batch1, batch2)

    countData_use <- limma::voom(in_data$countData, in_data$design1)
    fit <- limma::lmFit(countData_use$E, design = in_data$design1)
    fit <- limma::contrasts.fit(fit, in_data$contrast)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, number = nrow(countData_use$E))
    res$EventName <- rownames(res)

    res$AveExpr <- res$AveExpr - min(res$AveExpr)
    res <- as.data.table(res)
    
    rm(fit, in_data, countData_use)
    gc()
    return(res)
}

.ASE_edgeR_contrast <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
    in_data <- .ASE_contrast_expr(se, test_factor, 
        test_nom, test_denom,
        batch1, batch2)

    # countData_use <- limma::voom(in_data$countData, in_data$design1)
    y <- edgeR::DGEList(counts=in_data$countData, remove.zeros = FALSE)
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y,in_data$design1)
    
    fit <- edgeR::glmQLFit(y, in_data$design1)
    qlf <- edgeR::glmQLFTest(fit, contrast = in_data$contrast)

    res <- edgeR::topTags(qlf, n = nrow(y))
    res$table$EventName <- rownames(res)
    
    rm(fit, qlf, in_data, y)
    gc()
    return(as.data.table(res$table))
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
    
    fit <- edgeR::glmQLFit(y, model)
    return(list(
        fit = fit,
        model = model
    ))
}

.ASE_edgeR_fitASE <- function(se, model) {
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
    
    fit <- edgeR::glmQLFit(y, model)
    return(list(
        fit = fit,
        model = model
    ))
}

.ASE_limma_contrast_ts <- function(se, test_factor, # test_nom, test_denom,
        batch1, batch2, degrees) {
    in_data <- .ASE_contrast_expr_ts(
        se, test_factor, 
        # test_nom, test_denom,
        batch1, batch2, degrees
    )

    countData_use <- limma::voom(in_data$countData, in_data$design1)
    fit <- limma::lmFit(countData_use$E, design = in_data$design1)
    # fit <- limma::contrasts.fit(fit, in_data$contrast)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, coef = in_data$coef,
        number = nrow(countData_use$E))
    res$EventName <- rownames(res)

    res$AveExpr <- res$AveExpr - min(res$AveExpr)
    res <- as.data.table(res)
    
    rm(fit, in_data, countData_use)
    gc()
    return(res)
}


.ASE_limma_contrast_ASE <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
    in_data <- .ASE_contrast_ASE(se, test_factor, 
        test_nom, test_denom,
        batch1, batch2)

    countData_use <- limma::voom(in_data$countData, in_data$design1, 
        lib.size = 1)

    fit <- limma::lmFit(countData_use$E, design = in_data$design1)

    fit <- limma::contrasts.fit(fit, in_data$contrast)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, number = nrow(countData_use$E))
    res$EventName <- rownames(res)
    res$AveExpr <- res$AveExpr - min(res$AveExpr)
    res <- as.data.table(res)
    
    rm(fit, in_data, countData_use)
    gc()
    return(res)
}

.ASE_edgeR_contrast_ASE <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
    in_data <- .ASE_contrast_ASE(se, test_factor, 
        test_nom, test_denom,
        batch1, batch2)

    y <- edgeR::DGEList(counts=in_data$countData)
    y <- edgeR::estimateDisp(y,in_data$design1)
    
    fit <- edgeR::glmQLFit(y, in_data$design1)
    qlf <- edgeR::glmQLFTest(fit, contrast = in_data$contrast)

    res <- edgeR::topTags(qlf, n = nrow(y))
    res$table$EventName <- rownames(res)
    
    rm(fit, qlf, in_data, y)
    gc()
    return(as.data.table(res$table))
}

.ASE_limma_contrast_ASE_ts <- function(se, test_factor, # test_nom, test_denom,
        batch1, batch2, degrees) {
    in_data <- .ASE_contrast_ASE_ts(se, test_factor, 
        # test_nom, test_denom,
        batch1, batch2, degrees)

    countData_use <- limma::voom(in_data$countData, in_data$design1, 
        lib.size = 1)

    fit <- limma::lmFit(countData_use$E, design = in_data$design1)

    fit <- limma::contrasts.fit(fit, in_data$contrast)
    fit <- limma::eBayes(fit)

    res <- limma::topTable(fit, number = nrow(countData_use$E))
    res$EventName <- rownames(res)
    res$AveExpr <- res$AveExpr - min(res$AveExpr)
    res <- as.data.table(res)
    
    rm(fit, in_data, countData_use)
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
        batch1, batch2, BPPARAM) {
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
    if(BPPARAM$workers > 1) {
        sumExp <- satuRn::fitDTU(
            object = sumExp,
            formula = formula1,
            verbose = FALSE,
            parallel = TRUE,
            BPPARAM = BPPARAM
        )    
    } else {
        sumExp <- satuRn::fitDTU(
            object = sumExp,
            formula = formula1,
            verbose = FALSE,
            parallel = FALSE
        )
    }
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

# Adds human-readable labels
.ASE_add_flags <- function(res) {
    res_IR <- res[get("EventType") == "IR"]
    res_nonIR <- res[get("EventType") != "IR"]

    flags_IR <- rep("", nrow(res_IR))
    flags_nonIR <- rep("", nrow(res_nonIR))
    
    if(length(flags_IR) > 0) {
        flags_IR[res_IR$NMD_direction == 1] <- paste0(
            flags_IR[res_IR$NMD_direction == 1], ";Inc-NMD")
        flags_IR[res_IR$NMD_direction == -1] <- paste0(
            flags_IR[res_IR$NMD_direction == -1], ";Exc-NMD")
    }    
    if(length(flags_nonIR) > 0) {
        flags_nonIR[res_nonIR$NMD_direction == 1] <- paste0(
            flags_nonIR[res_nonIR$NMD_direction == 1], ";Inc-NMD")
        flags_nonIR[res_nonIR$NMD_direction == -1] <- paste0(
            flags_nonIR[res_nonIR$NMD_direction == -1], ";Exc-NMD")
            
        isIncNovel <- grepl("novel", tstrsplit(res_nonIR$EventName, 
            split = ";", fixed = TRUE)[[1]])
        isExcNovel <- grepl("novel", tstrsplit(res_nonIR$EventName, 
            split = ";", fixed = TRUE)[[2]])
        flags_nonIR[isIncNovel] <- paste0(
            flags_nonIR[isIncNovel], ";Inc-novel")
        flags_nonIR[isExcNovel] <- paste0(
            flags_nonIR[isExcNovel], ";Exc-novel")
    }
    flagsStr <- c(flags_IR, flags_nonIR)
    res[, c("flags") := substr(flagsStr, 2, nchar(flagsStr))]
    return(res[, c("EventName","EventType","EventRegion", "flags")])
}

.ASE_add_diag <- function(res, se, test_factor, test_nom, test_denom) {
    rowData <- as.data.frame(rowData(se))
    rowData.DT <- .ASE_add_flags(as.data.table(rowData[,
        c("EventName","EventType","EventRegion", "NMD_direction")]))
        
    diag <- makeMeanPSI(se, res$EventName,
        test_factor, list(test_nom, test_denom))
    colnames(diag)[2:3] <- c(paste0("AvgPSI_", test_nom),
        paste0("AvgPSI_", test_denom))
    diag$deltaPSI <- diag[, 2] - diag[, 3]
    if(all(c("EventName","EventType","EventRegion", "flags") %in%
        colnames(res))) {
        res <- cbind(
            res[,c("EventName")],
            as.data.table(round(diag[,-1], 4)),
            res[,-c("EventName","EventType","EventRegion", "flags")]
        )            
    } else {
        res <- cbind(
            res[,c("EventName")],
            as.data.table(round(diag[,-1], 4)),
            res[,-c("EventName")]
        )
    }
    res <- rowData.DT[res, on = "EventName"]
    return(res)
}

.ASE_add_diag_multi <- function(
        res, se, test_factor, 
        conditionList
) {
    rowData <- as.data.frame(rowData(se))
    rowData.DT <- .ASE_add_flags(as.data.table(rowData[,
        c("EventName","EventType","EventRegion", "NMD_direction")]))
    diag <- makeMeanPSI(se, res$EventName,
        test_factor, conditionList)
    for(i in seq_len(length(conditionList))) {
        colnames(diag)[i+1] <- paste0("AvgPSI_", conditionList[i])
    }
    if(ncol(diag) == 3) {
        diag$deltaPSI <- diag[, 2] - diag[, 3]
    }
    if(all(c("EventName","EventType","EventRegion", "flags") %in%
        colnames(res))) {
        res <- cbind(
            res[,c("EventName")],
            as.data.table(round(diag[,-1], 4)),
            res[,-c("EventName","EventType","EventRegion", "flags")]
        )            
    } else {
        res <- cbind(
            res[,c("EventName")],
            as.data.table(round(diag[,-1], 4)),
            res[,-c("EventName")]
        )
    }
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

################################################################################

.ASE_contrast_expr <- function(
        se, test_factor, 
        test_nom, test_denom,
        batch1, batch2
) {
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
    
    return(list(
        design1 = design1,
        countData = countData,
        contrast = contrast
    ))
}

.ASE_contrast_expr_ts <- function(
        se, test_factor, 
        # test_nom, test_denom,
        batch1, batch2,
        degrees = 3
) {
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
    if(degrees > length(unique(colData[, test_factor])) - 1)
        degrees <- length(unique(colData[, test_factor])) - 1
    if(degrees < 1) {
        .log(paste(
            "Not enough degrees of freedom in", test_factor, "variable.",
            "Perhaps just use a simple contrast between conditions?"
        ))
    }

    condition_factor <- splines::ns(colData[, test_factor], 
        df = degrees)
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
    coef <- seq(ncol(design1) - degrees + 1, ncol(design1))

    return(list(
        design1 = design1,
        countData = countData,
        coef = coef
    ))
}

.ASE_contrast_ASE <- function(
        se, test_factor, 
        test_nom, test_denom,
        batch1, batch2
) {
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
    
    return(list(
        design1 = design1,
        countData = countData,
        contrast = contrast
    ))
}

.ASE_contrast_ASE_ts <- function(
        se, test_factor, 
        batch1, batch2,
        degrees = 3
) {
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

    if(degrees > length(unique(colData[, test_factor])) - 1)
        degrees <- length(unique(colData[, test_factor])) - 1
    if(degrees < 1) {
        .log(paste(
            "Not enough degrees of freedom in", test_factor, "variable.",
            "Perhaps just use a simple contrast between conditions?"
        ))
    }
    condition_factor <- splines::ns(colData[, test_factor], 
        df = degrees)
    ASE <- colData[, "ASE"]
    if(batch2 != "") {
        batch2_factor <- colData[, batch2]
        batch1_factor <- colData[, batch1]
        design1 <- model.matrix(~0 + batch1_factor + batch2_factor + 
            ASE + ASE:condition_factor)
    } else if(batch1 != "") {
        batch1_factor <- colData[, batch1]
        design1 <- model.matrix(~0 + batch1_factor + 
            ASE + ASE:condition_factor)
    } else {
        design1 <- model.matrix(~0 + ASE + ASE:condition_factor)
    }
    colnames(design1) <- sub(":",".",colnames(design1))
    contrast <- rep(0, ncol(design1))
    contrast[which(grepl("ASEExcluded.", colnames(design1)))] <- -1
    contrast[which(grepl("ASEIncluded.", colnames(design1)))] <- 1
    
    return(list(
        design1 = design1,
        countData = countData,
        contrast = contrast
    ))
}