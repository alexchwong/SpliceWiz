# NxtSE Class functions

#' The NxtSE class
#'
#' The NxtSE class inherits from the \linkS4class{SummarizedExperiment} 
#' class and is constructed from [makeSE]. NxtSE extends SummarizedExperiment
#' by housing additional assays pertaining to IR and splice junction counts.
#' @param x A NxtSE object
#' @param i,j Row and column subscripts to subset a NxtSE object.
#' @param ... In NxtSE(), additional arguments to be passed onto
#'   SummarizedExperiment()
#' @param withDimnames (default TRUE) Whether exported assays should be
#'   supplied with row and column names of the NxtSE object.
#'   See \linkS4class{SummarizedExperiment}
#' @param deparse.level See [base::cbind] for a description of this argument.
#' @param drop A logical(1), ignored by these methods.
#' @param value The value to replace. Must be a matrix for the 
#'   up_inc<-, down_inc<-, up_exc<- and down_exc<- replacers, 
#'   and a character vector for covfile<-
#' @return See Functions section (below) for details
#' @examples
#'
#' # Run the full pipeline to generate a NxtSE object:
#'
#' buildRef(
#'     reference_path = file.path(tempdir(), "Reference"),
#'     fasta = chrZ_genome(), 
#'     gtf = chrZ_gtf()
#' )
#'
#' bams <- SpliceWiz_example_bams()
#' processBAM(bams$path, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "SpliceWiz_Output")
#' )
#' 
#' expr <- findSpliceWizOutput(file.path(tempdir(), "SpliceWiz_Output"))
#' collateData(expr, 
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "Collated_output")
#' )
#' 
#' se <- makeSE(collate_path = file.path(tempdir(), "Collated_output"))
#'
#' # Coerce NxtSE -> SummarizedExperiment
#' se_raw <- as(se, "SummarizedExperiment")
#' 
#' # Coerce SummarizedExperiment -> NxtSE
#' se_NxtSE <- as(se_raw, "NxtSE")
#' identical(se, se_NxtSE) # Returns TRUE
#' 
#' # Get Junction reads of SE / MXE and spans-reads of IR events
#' up_inc(se)
#' down_inc(se)
#' up_exc(se)
#' down_exc(se)
#' 
#' # Get list of available coverage files
#' covfile(se)
#' 
#' # Get sample QC information
#' sampleQC(se)
#'
#' # Get resource NxtIRF data (used internally for plotCoverage())
#' cov_data <- ref(se)
#' names(cov_data)
#'
#' # Subset functions
#' se_by_samples <- se[,1:3]
#' se_by_events <- se[1:10,]
#' se_by_rowData <- subset(se, EventType == "IR")
#'
#' # Cbind (bind event_identical NxtSE by samples)
#' se_by_samples_1 <- se[,1:3]
#' se_by_samples_2 <- se[,4:6]
#' se_cbind <- cbind(se_by_samples_1, se_by_samples_2)
#' identical(se, se_cbind) # should return TRUE
#'
#' # Rbind (bind sample_identical NxtSE by events)
#' se_IR <- subset(se, EventType == "IR")
#' se_SE <- subset(se, EventType == "SE")
#' se_IRSE <- rbind(se_IR, se_SE)
#' identical(se_IRSE, subset(se, EventType %in% c("IR", "SE"))) # TRUE
#'
#' # Convert HDF5-based NxtSE to in-memory se
#' # makeSE() creates a HDF5-based NxtSE object where all assay data is stored
#' # as an h5 file instead of in-memory. All operations are performed as
#' # delayed operations as per DelayedArray package.
#' # To realize the NxtSE object as an in-memory object, use:
#' 
#' se_real <- realize_NxtSE(se)
#' identical(se, se_real) # should return FALSE
#'
#' # To check the difference, run:
#' class(up_inc(se))
#' class(up_inc(se_real))
#'
#' @name NxtSE-class
#' @aliases
#' NxtSE-methods
#' up_inc up_inc,NxtSE-method
#' up_inc<- up_inc<-,NxtSE-method
#' down_inc down_inc,NxtSE-method
#' down_inc<- down_inc<-,NxtSE-method
#' up_exc up_exc,NxtSE-method
#' up_exc<- up_exc<-,NxtSE-method
#' down_exc down_exc,NxtSE-method
#' down_exc<- down_exc<-,NxtSE-method
#' covfile covfile,NxtSE-method
#' covfile<- covfile<-,NxtSE-method
#' sampleQC sampleQC,NxtSE-method
#' sampleQC<- sampleQC<-,NxtSE-method
#' ref ref,NxtSE-method
#' realize_NxtSE realize_NxtSE,NxtSE-method
#' coerce,SummarizedExperiment,NxtSE-method
#' @md
#' @export
setClass("NxtSE",
    slots=c(int_elementMetadata = "DataFrame",
        int_colData = "DataFrame",
        int_metadata = "list"),
    contains = "SummarizedExperiment"
)


#' NxtIRF filters to remove low-abundance alternative splicing and intron
#' retention events
#'
#' @param filterClass Must be either `"Data"` or `"Annotation"`. See details
#' @param filterType If `filterClass = "Data"`, then must be one of 
#'   `c("Depth", "Coverage", "Consistency")`. If `filterClass = "Annotation"`,
#'   must be one of `c("Protein_Coding", "NMD", "TSL")`. See details
#' @param pcTRUE If conditions are set, what percentage of all samples in each
#'   of the condition must the filter be satisfied for the event to pass the
#'   filter check. Must be between 0 and 100 (default 100)
#' @param minimum Filter-dependent argument. See details
#' @param maximum Filter-dependent argument. See details
#' @param minDepth Filter-dependent argument. See details
#' @param condition (default "") If set, must match the name of an experimental
#'   condition in the NxtSE object to be filtered, 
#'   i.e. a column name in `colData(se)`. Leave blank to disable filtering
#'   by condition
#' @param minCond (default -1) If condition is set, how many minimum number of
#'   conditions must pass the filter criteria. For example,
#'   if condition = "Batch", and batches are "A", "B", or "C", setting
#'   `minCond = 2` with `pcTRUE = 100` means that all samples belonging to
#'   two of the three types of `Batch` must pass the filter criteria.
#'   Setting `-1` means all elements of `condition` must
#'   pass criteria. Set to `-1` when the number of elements in the experimental
#'   condition is unknown. Ignored if `condition` is left blank.
#' @param EventTypes What types of events are considered for filtering. Must be 
#'   one of `c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")`. Events 
#'   not specified in `EventTypes` are not filtered (i.e. they will pass the
#'   filter without checks)
#' @details
#'   **Annotation Filters**
#'   * **Protein_Coding**: Filters for alternative splicing or IR events 
#'       involving protein-coding transcripts.
#'       No additional parameters required.
#'   * **NMD**: Filters for events in which one isoform is a
#'       predicted NMD substrate.
#'   * **TSL**: filters for events in which both
#'       isoforms have a TSL level below or equal to `minimum`
#'   * **Terminus** (New as of version 1.1.1): 
#'       In alternate first exons, the splice junction must
#'       not be shared with another transcript for which it is not its first
#'       intron. For alternative last exons, the splice junction must not be
#'       shared with another transcript for which it is not its last intron
#'   * **ExclusiveMXE** (New as of version 1.1.1):
#'       For MXE events, the two alternate
#'       casette exons must not overlap in their genomic regions
#'
#'   **Data Filters**
#'   * **Depth**: Filters IR or alternative splicing events of transcripts
#'       that are "expressed" with adequate `Depth` as calculated by the
#'       sum of all splicing and IR reads spanning the event. Events with 
#'       `Depth` below `minimum` are filtered out
#'   * **Coverage**: Coverage means different things to IR and alternative
#'       splicing.\cr\cr
#'     For **IR**, Coverage refers to the percentage of the measured intron
#'       covered with reads. Introns of samples with an IntronDepth above 
#'       `minDepth` are assessed, with introns with coverage 
#'       below `minimum` are filtered out.\cr\cr
#'     For **Alternative Splicing**, Coverage refers to the percentage of all
#'       splicing events observed across the genomic region that is compatible
#'       with either the included or excluded event. This prevents NxtIRF from 
#'       doing differential analysis between two minor isoforms. Instead of
#'       IntronDepth, in AS events NxtIRF considers events where the spliced
#'       reads from both exonic regions exceed `minDepth`.
#'       Then, events with a splicing coverage below `minimum`
#'       are excluded. \cr\cr
#'       We recommend testing IR events for > 90% coverage and AS
#'       events for > 60% coverage as given in the default filters which can be
#'       accessed using [getDefaultFilters]\cr\cr
#'    * **Consistency**: Skipped exons (SE) and mutually exclusive exons
#'       (MXE) comprise reads aligned to two contiguous splice junctions. 
#'       Most algorithms take the average counts from both junctions. This
#'       will inadvertently include transcripts that share one but not both
#'       splice events. To check that this is not happening, we require both
#'       splice junctions to have comparable counts.
#'       This filter checks whether reads from each splice junction comprises
#'       a reasonable proportion of the sum of these reads.\cr\cr
#'     Events are excluded if either of the upstream or downstream
#'       event is lower than total splicing events by a log-2 magnitude 
#'       above `maximum`. For example, if 
#'       `maximum = 2`, we require both upstream and downstream
#'       events to represent at least 1/(2^2) = 1/4 of the sum of upstream
#'       and downstream event. If `maximum = 3`, then each junction must be at
#'       least 1/8 of total, etc.
#'       This is considered for each isoform of each event, as long as the
#'       total counts belonging to the considered isoform is above 
#'       `minDepth`.\cr\cr
#'     IR-events are also checked. For IR events, the upstream and downstream
#'     exon-intron spanning reads must comprise a reasonable proportion of total
#'     exon-intron spanning reads.
#'
#'   We highly recommend using the default filters, which can be acquired 
#'     using [getDefaultFilters]
#' @return A ASEFilter object with the specified parameters
#' @examples
#' # Create a ASEFilter that filters for protein-coding ASE
#' f1 <- ASEFilter(filterClass = "Annotation", filterType = "Protein_Coding")
#'
#' # Create a ASEFilter that filters for Depth >= 20 in IR events
#' f2 <- ASEFilter(
#'     filterClass = "Data", filterType = "Depth",
#'     minimum = 20, EventTypes = c("IR", "RI")
#' )
#'
#' # Create a ASEFilter that filters for Coverage > 60% in splice events
#' # that must be satisfied in at least 2 categories of condition "Genotype"
#' f3 <- ASEFilter(
#'     filterClass = "Data", filterType = "Coverage",
#'     minimum = 60, EventTypes = c("MXE", "SE", "AFE", "ALE", "A3SS", "A5SS"),
#'     condition = "Genotype", minCond = 2
#' )
#'
#' # Create a ASEFilter that filters for Depth > 10 in all events
#' # that must be satisfied in at least 50% of each gender
#' f4 <- ASEFilter(
#'     filterClass = "Data", filterType = "Depth",
#'     minimum = 10, condition = "gender", pcTRUE = 50
#' )
#'
#' # Get a description of what these filters do:
#' f1
#' f2
#' f3
#' f4
#'
#' @name ASEFilter-class
#' @aliases ASEFilter
#' @seealso [Run_SpliceWiz_Filters]
#' @md
#' @export
setClass("ASEFilter",
    slots = c(
        filterClass = "character",
        filterType = "character",
    
        pcTRUE = "numeric",
        minimum = "numeric",
        maximum = "numeric",
        minDepth = "numeric",
        condition = "character",
        minCond = "numeric",
        EventTypes = "character"
    )
)