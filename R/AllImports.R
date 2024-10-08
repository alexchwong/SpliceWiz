#' @useDynLib SpliceWiz, .registration = TRUE
#' @import NxtIRFdata
#' @import ompBAM
#' @importFrom methods as is coerce callNextMethod new validObject show
#' @importFrom parallel detectCores
#' @importFrom scales label_number cut_si
#' @importFrom stats as.formula model.matrix qt runif na.omit prcomp aggregate
#' @importFrom utils download.file packageVersion getFromNamespace
#' @importFrom utils setTxtProgressBar txtProgressBar head
#' @importFrom tools R_user_dir
#' @importFrom magrittr %>%
#' @importFrom R.utils gzip
#' @importFrom Rcpp evalCpp
#' @import data.table
#' @importFrom fst read.fst write.fst threads_fst
#' @import ggplot2
#' @importFrom AnnotationHub AnnotationHub cache query removeResources
#' @importFrom RSQLite SQLite
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcquery
#' @importFrom BiocGenerics start end start<- end<- width
#' @importFrom BiocGenerics nrow ncol rbind cbind
#' @importFrom BiocParallel SnowParam MulticoreParam SerialParam
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom Biostrings readDNAStringSet DNAStringSet 
#' @importFrom Biostrings translate
#' @importFrom stringi stri_locate_first_fixed
#' @importFrom BSgenome getSeq
#' @importFrom DelayedArray qlogis plogis rowMeans DelayedArray rowSums
#' @importFrom DelayedMatrixStats rowSds rowAnyNAs rowMaxs
#' @importFrom genefilter rowttests
#' @importFrom GenomeInfoDb sortSeqlevels seqinfo seqlengths seqlevels<- 
#' @importFrom GenomeInfoDb seqlevels seqlevelsStyle seqlevelsStyle<-
#' @importFrom GenomeInfoDb genomeStyles seqinfo<-
#' @importFrom GenomicRanges GRanges reduce findOverlaps 
#' @importFrom GenomicRanges makeGRangesFromDataFrame 
#' @importFrom GenomicRanges makeGRangesListFromDataFrame split strand 
#' @importFrom GenomicRanges flank setdiff seqnames psetdiff disjoin
#' @importFrom GenomicRanges strand<- seqnames<-
#' @importFrom HDF5Array HDF5Array writeHDF5Array h5writeDimnames 
#' @importFrom IRanges IRanges Views RleList
#' @importFrom progress progress_bar
#' @importFrom patchwork area plot_layout wrap_plots
#' @importFrom plotly config layout plotlyOutput event_data ggplotly 
#' @importFrom plotly plotlyProxy plotlyProxyInvoke renderPlotly subplot 
#' @importFrom plotly highlight orca toWebGL partial_bundle
#' @importFrom plotly plot_ly add_trace add_ribbons plotly_build
#' @importFrom plotly plot_ly event_register event_unregister
#' @importFrom rhdf5 h5createFile h5createDataset h5delete h5write h5createGroup
#' @importFrom rtracklayer import export TwoBitFile track
#' @importFrom S4Vectors mcols mcols<- metadata Rle metadata<- SimpleList 
#' @importFrom S4Vectors from to setValidity2 DataFrame
#' @importFrom S4Vectors bindCOLS bindROWS getListElement setListElement
#' @importFrom SummarizedExperiment SummarizedExperiment 
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<-
#' @importFrom SummarizedExperiment assay assays assay<- assays<-
#' @importFrom SummarizedExperiment assayNames assayNames<-
#' @importClassesFrom S4Vectors DataFrame 
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import shiny
#' @import shinydashboard
#' @importFrom shinyFiles getVolumes parseDirPath parseFilePaths parseSavePath
#' @importFrom shinyFiles shinyDirButton shinyDirChoose shinyFileChoose
#' @importFrom shinyFiles shinyFilesButton shinyFileSave shinySaveButton
#' @importFrom shinyWidgets sliderTextInput updateSliderTextInput
#' @importFrom shinyWidgets prettyRadioButtons materialSwitch
#' @importFrom shinyWidgets radioGroupButtons updateRadioGroupButtons
#' @importFrom shinyWidgets switchInput actionBttn updateSwitchInput
#' @importFrom shinyWidgets sendSweetAlert ask_confirmation
#' @importFrom rhandsontable rhandsontable hot_to_r hot_col
#' @importFrom rhandsontable renderRHandsontable rHandsontableOutput
#' @importFrom htmltools tagList singleton
#' @importFrom stats runif na.omit prcomp complete.cases p.adjust
#' @importFrom DT datatable selectRows dataTableProxy  
#' @importFrom grDevices colorRampPalette
#' @importFrom heatmaply heatmaply
#' @importFrom pheatmap pheatmap
#' @importFrom matrixStats colVars
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom rvest read_html html_nodes html_attr
#' @importFrom httr GET timeout
NULL

# Import namespaces of internal functions of BiocGenerics and S4Vectors
# Used to set up NxtSE class

BG_replaceSlots <- getFromNamespace("replaceSlots", "BiocGenerics")
S4_disableValidity <- getFromNamespace("disableValidity", "S4Vectors")
S4_selectSome <- getFromNamespace("selectSome", "S4Vectors")

sWidgets_attachDep <- getFromNamespace("attachShinyWidgetsDep", "shinyWidgets")
sFiles_formatFiletype <- getFromNamespace("formatFiletype", "shinyFiles")

# Checks character indices on NxtSE object
SE_charbound <- function(idx, txt, fmt) {
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
        msg <- paste(S4_selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}
