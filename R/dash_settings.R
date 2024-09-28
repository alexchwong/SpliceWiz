# 
setreactive_system <- function() {
    reactiveValues(
        threads_initialized = FALSE,
        n_threads = 1
    )
}

# Settings for "Make New SpliceWiz Reference"
setreactive_newref <- function() {
    reactiveValues(
        newref_path = "",
        newref_fasta = "",
        newref_gtf = "",
        newref_mappa = "",
        newref_NPA = "",
        newref_bl = "",

    # FTP data
        ensemblWorking = NULL,
        availRelease = "",
        availSpecies = "",
        availGenome = "",
        availGTF = ""
    )
}

# Settings for "Make New SpliceWiz Experiment"
setreactive_expr <- function() {
    reactiveValues(
        ref_path = "",
        bam_path = "",
        sw_path = "",
        anno_file = "",
        collate_path = "",
        collate_path_prompt = "",

        df.bams = c(),
        df.files = c(),
        df.anno = c(),
        disallow_df_update = FALSE,

        df.files_savestate = c(),
        df.anno_savestate = c(),

        ref_settings = c(),
        ref_table = data.frame(),
        
        expr_path = "",
        selected_bams = c(),
        collateData_args = c(),
        df = c(),
        se = NULL,
        
        se_expr_infobox_state = -1
    )
}

# Settings for navQC
setreactive_QC <- function() {
    reactiveValues(
        QC = c()
    )
}

# Settings for navFilter
setreactive_filtered_SE <- function() {
    reactiveValues(
        filterSummary = NULL,
        filters = list(),
        filtered_SE = NULL
    )
}

# Settings for navAnalyse
setreactive_DE <- function() {
    reactiveValues(
        res = NULL,
        res_settings = list(),
        method = NULL,
        batchVar1 = NULL,
        batchVar2 = NULL,
        DE_Var = NULL,
        nom_DE = NULL,
        denom_DE = NULL,
        IRmode_DE = NULL,
        dof = 1,
        
        DT_DE_rows_all = NULL,
        DT_DE_rows_selected = NULL,
        command_selected = NULL,
        filters = list()
    )
}

setreactive_GO <- function() {
    reactiveValues(
        trigger = NULL,
        
        res = NULL,

        errorMsg = "",

        # gene_ids that can optionally be saved to file by user
        gene_ids = NULL,
        univ_ids = NULL,

        # Save GO results
        filteredVolc = NULL,
        
        resGO = NULL,

        final_plot = NULL,
        ggplot = NULL
    )
}

# Settings for Diag and Volcano
setreactive_Diag <- function() {
    # NB same code as Volcano
    reactiveValues(
        trigger = NULL,
        
        useDE = NULL,
        
        meanPSI = NULL,

        plot_ini = FALSE,
        plotly_click = NULL,
        final_plot = NULL,
        ggplot = NULL,
        selected = NULL
    )
}

# Settings for Heatmap
setreactive_Heat <- function() {
    reactiveValues(
        trigger = NULL,
        
        useDE = NULL,
        eventsGO = NULL,
        mat = NULL,
        
        ggplot = NULL,
        final_plot = NULL
    )
}

# Settings for Coverage plots
setreactive_Cov <- function() {
    reactiveValues(
        geneList = NULL,
        useDE = NULL,
        
        view_chr = "",
        view_start = "",
        view_end = "",
        data_start = 0,
        data_end = 0,

        view_strand = "*",

        event.ranges = NULL,

        plotly_relayout = NULL,
        plot_ini = FALSE,
        
        final_plot = NULL,
        
        trigger = NULL,
        plot_params = NULL
    )
}

# Settings for Coverage plots (NEW)
setreactive_Cov2 <- function() {
    reactiveValues(
        # geneList = NULL,
        useDE = NULL,
        
        trackTable = data.frame(),
        exonsTable = data.frame(),
        transcripts = data.frame(),
        exons_gr = GRanges(),

    # New ranges can be triggered from different sources
        newGR = GRanges(), # aggregate GRanges
    # This is used to trigger plot refresh in absence of locale change
        plotTrigger = NULL,
        
        dataObj = covDataObject(),
        plotObj = covPlotObject(),
        plotlyObj = covPlotly(),
        plotlyFig = plot_ly(),
        plotCount = 0,

        event.ranges = NULL,
        prevEventGR = NULL,

        plotly_relayout = NULL,
        # plot_ini = FALSE,
        oldPlotSettings = list(),
        prevReqEvent = NULL,
        plotReq = NULL,
        normEvent = list(),
        
        ggplot = ggplot(),
        exonsplot = ggplot()
    )
}