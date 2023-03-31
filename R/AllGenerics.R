# Generics from SummarizedExperiment

#' @export
setGeneric("rowData", signature="x",
    function(x, use.names=TRUE, ...) standardGeneric("rowData"))

#' @export
setGeneric("rowData<-",
    function(x, ..., value) standardGeneric("rowData<-"))
    
#' @export
setGeneric("colData", function(x, ...) standardGeneric("colData"))

#' @export
setGeneric("colData<-",
    function(x, ..., value) standardGeneric("colData<-"))
    
#' @export
setGeneric("assays", signature="x",
    function(x, withDimnames=TRUE, ...) standardGeneric("assays"))

#' @export
setGeneric("assays<-", signature=c("x", "value"),
    function(x, withDimnames=TRUE, ..., value) standardGeneric("assays<-"))

#' @export
setGeneric("assay", signature=c("x", "i"),
    function(x, i, withDimnames=TRUE, ...) standardGeneric("assay"))

#' @export
setGeneric("assay<-", signature=c("x", "i"),
    function(x, i, withDimnames=TRUE, ..., value) standardGeneric("assay<-"))

#' @export
setGeneric("assayNames", function(x, ...) standardGeneric("assayNames"))

#' @export
setGeneric("assayNames<-",
    function(x, ..., value) standardGeneric("assayNames<-"))

# NxtSE specific functions:

setGeneric("realize_NxtSE", 
    function(x, includeJunctions = FALSE, withDimnames=TRUE, ...)
    standardGeneric("realize_NxtSE"))

setGeneric("update_NxtSE", function(x, ...) standardGeneric("update_NxtSE"))

setGeneric("up_inc", 
    function(x, withDimnames=TRUE, ...) standardGeneric("up_inc"))

setGeneric("down_inc", 
    function(x, withDimnames=TRUE, ...) standardGeneric("down_inc"))

setGeneric("up_exc", 
    function(x, withDimnames=TRUE, ...) standardGeneric("up_exc"))

setGeneric("down_exc", 
    function(x, withDimnames=TRUE, ...) standardGeneric("down_exc"))

setGeneric("covfile", 
    function(x, withDimnames=TRUE, ...) standardGeneric("covfile"))

setGeneric("sampleQC", 
    function(x, withDimnames=TRUE, ...) standardGeneric("sampleQC"))

setGeneric("sourcePath", 
    function(x, withDimnames=TRUE, ...) standardGeneric("sourcePath"))

setGeneric("ref", 
    function(x, withDimnames=TRUE, ...) standardGeneric("ref"))

setGeneric("junc_PSI", 
    function(x, withDimnames=TRUE, ...) standardGeneric("junc_PSI"))

setGeneric("junc_counts", 
    function(x, withDimnames=TRUE, ...) standardGeneric("junc_counts"))

setGeneric("junc_counts_uns", 
    function(x, withDimnames=TRUE, ...) standardGeneric("junc_counts_uns"))

setGeneric("row_gr", 
    function(x, withDimnames=TRUE, ...) standardGeneric("row_gr"))

setGeneric("junc_gr", 
    function(x, withDimnames=TRUE, ...) standardGeneric("junc_gr"))

setGeneric("up_inc<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("up_inc<-"))

setGeneric("down_inc<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("down_inc<-"))

setGeneric("up_exc<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("up_exc<-"))

setGeneric("down_exc<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("down_exc<-"))

setGeneric("covfile<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("covfile<-"))

setGeneric("sampleQC<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("sampleQC<-"))
    
setGeneric("ref<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("ref<-"))

setGeneric("sourcePath<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("sourcePath<-"))

setGeneric("row_gr<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("row_gr<-"))
    
setGeneric("junc_PSI<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("junc_PSI<-"))
    
setGeneric("junc_counts<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("junc_counts<-"))

setGeneric("junc_counts_uns<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("junc_counts_uns<-"))
    
setGeneric("junc_gr<-",
    function(x, withDimnames=TRUE, ..., value) standardGeneric("junc_gr<-"))
    
    
# covPlotObject specific functions:

setGeneric("tracks", signature="x", function(
    x, 
    ...
) standardGeneric("tracks"))

setGeneric("condition", signature="x", function(
    x, 
    ...
) standardGeneric("condition"))

setGeneric("getExonRanges", signature="x", function(
    x, 
    ...
) standardGeneric("getExonRanges"))
